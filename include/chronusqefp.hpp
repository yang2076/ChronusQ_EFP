#ifndef __INCLUDED_CHRONUSQ_HEADER_HPP__
#define __INCLUDED_CHRONUSQ_HEADER_HPP__

#include <util/matout.hpp>
#include <string>
#include <memmanager.hpp>
#include "private.h"
#include <aointegrals.hpp>
#include <singleslater.hpp>
#include <cxxapi/input.hpp>
namespace ChronusQ{

  class EFPBase{
  public:
    int type_choice;
    EFPBase() = default;
    EFPBase(const EFPBase &) = default;
    EFPBase(EFPBase &&)      = default;
    
    virtual void Initialize(::efp_opts* option, 
                           std::string potential_file, 
                           struct Fragment_ifm* frag_ifm) = 0;

    virtual void EFP_Compute(int do_gradient) = 0;
    virtual void Wavefunction_dependent(double* energy_ptr) = 0; 

  };
  
  template<typename IntsT, typename MatsT> 
  class EFP : public EFPBase{
  private:
    static ::efp_result elec_field_fn(size_t n_pt, const double* pt_coord, double* field, void* user_data){
       auto ss_1 = static_cast<SingleSlaterBase*>(user_data);
       auto slater0 = dynamic_cast<SingleSlater<double,double>* >(ss_1);
       auto slater1 = dynamic_cast<SingleSlater<double,dcomplex>* >(ss_1);
       auto slater2 = dynamic_cast<SingleSlater<dcomplex,double>* >(ss_1);
       auto slater3 = dynamic_cast<SingleSlater<dcomplex,dcomplex>* >(ss_1);

       for(int iPol = 0; iPol < n_pt; iPol++){
         for(int iXYZ = 0; iXYZ < 3; iXYZ++){
           if(slater0 != NULL){
           field[iPol*3+iXYZ] = -1*slater0->template computeOBProperty<double, SCALAR>(slater0->aoints.indE1[iPol][iXYZ]);
           }
           if(slater1 != NULL){
           field[iPol*3+iXYZ] = -1*slater1->template computeOBProperty<double, SCALAR>(slater1->aoints.indE1[iPol][iXYZ]);
           }
           if(slater2 != NULL){
           field[iPol*3+iXYZ] = -1*slater2->template computeOBProperty<double, SCALAR>(slater2->aoints.indE1[iPol][iXYZ]);
           }
           if(slater3 != NULL){
           field[iPol*3+iXYZ] = -1*slater3->template computeOBProperty<double, SCALAR>(slater3->aoints.indE1[iPol][iXYZ]);
           }
         }
       }
      return EFP_RESULT_SUCCESS;
    };  
    static size_t          n_pt;
    static const double*   pt_coord;
    static double*         field;

  protected:    
    CQMemManager    &memManager_;
    enum efp_result result_;

    
    enum efp_coord_type coord_type; 
    
    size_t          n_pc;
    const double*   pc;
    const double*   pc_coord;
    const double*   frag_coord;
    size_t          frag_idx;
    double*         charge;
    size_t          n_atom;
    size_t          n_frag;
    size_t          n_mult;
    size_t 	    n_frag_mult;
    double*         mult_coord;
    double*         mult;
    size_t          n_pol;
    double*         pol_coord;
    double*         pol;
    double*         pol_conj;
    double*         energy_denp;
    double*         grad;
    double*         grad_pc;
    AOIntegrals<IntsT> &aoints;
    SingleSlaterBase*  ss_;
    SingleSlaterBase*  propagator_;
  public:
    ::efp*          efp_;
    ::efp_energy    efp_energy_;
    std::vector<size_t> atom_num;
    std::vector<size_t> frag_mult_num;
    std::vector<double> atom_znuc;
    //ctors(Default copy and move)    
    EFP( const EFP & ) = default;
    EFP( EFP && )      = default;
    // No default ctor
    EFP() = delete;
    // Ctor needs integrals, single slater, and memManager
    EFP(AOIntegrals<IntsT>& ao, SingleSlaterBase* ss, CQMemManager& memManager):
		   aoints(ao), ss_(ss), memManager_(memManager), EFPBase(){};
    
       
    // Initialization
    void Initialize(::efp_opts* option, 
                    std::string potential_file, 
                    struct Fragment_ifm* frag_ifm){
      
      // print the banner of EFP 
      efp_print_banner();
      
      // create the EFP struct
      efp_ = efp_create();
      
      // Set the options and EFP file
      result_ = efp_set_opts(efp_, option);
      result_ = efp_add_potential(efp_, potential_file.c_str());
      coord_type = frag_ifm->Coord_type;
      frag_coord = frag_ifm->Frag_coord;
      for(std::vector<std::string>::iterator frag = (frag_ifm->fragments).begin();
          frag != (frag_ifm->fragments).end(); ++frag){
         result_ = efp_add_fragment(efp_, frag->c_str());
      }
      
      // Obtain the EFP fragment nuclei information
      result_ = efp_get_frag_count(efp_, &n_frag);
      n_atom = 0;
      for(int i = 0; i < n_frag; i++){
        size_t trans_;
        result_ = efp_get_frag_atom_count(efp_,i,&trans_);
        atom_num.push_back(trans_);
        n_atom += atom_num[i];
        result_ = efp_get_frag_multipole_count(efp_,i,&trans_);
        frag_mult_num.push_back(trans_);
      }
      int t = 0;
      for(int i = 0; i < n_frag; i++){
        auto efp_atom_trans = memManager_.template malloc<::efp_atom>(atom_num[i]);
        result_ = efp_get_frag_atoms(efp_,i,atom_num[i],efp_atom_trans);
        for(int j = 0; j < atom_num[i]; j++){
          atom_znuc.push_back((efp_atom_trans+j)->znuc);
          t++; 
        }
      }
      result_ = efp_set_coordinates(efp_, coord_type, frag_coord);

      // Obtain the ab initio nuclei information
      set_qm_nuclei();

      // Set the user data and field passing function
      efp_->get_electron_density_field_user_data = ss_;
      efp_->get_electron_density_field = &elec_field_fn;

      // Begin collecting the polarization information
      efp_get_induced_dipole_count(efp_,&n_pol);
      pol_coord = memManager_.template malloc<double>(3*n_pol);
      efp_get_induced_dipole_coordinates(efp_,pol_coord);
      
      // Begin collecting the multipole information
      efp_get_multipole_count(efp_, &n_mult);
      mult_coord = memManager_.template malloc<double>(3*n_mult);
      mult = memManager_.template malloc<double>(20*n_mult);
      
      efp_get_multipole_coordinates(efp_, mult_coord);
      efp_get_multipole_values(efp_, mult);
      int k = 0;
      t = 0;
      for(int i = 0; i < n_frag; i++){
        for(int j = 0; j < atom_num[i]; j++){
          *(mult+t) += atom_znuc[k];
          k++;
          t += 20;
        }
        t += (frag_mult_num[i]-atom_num[i])*20;
      }

      efp_prepare(efp_);
      aoints.computeEFP_integrals(pol_coord,&n_pol,1);

    };



    // Convey the ab initio nuclei charge to EFP system.
    void set_qm_nuclei(){
      auto n_pc = aoints.molecule_.nAtoms;
      auto pc_ = memManager_.template malloc<double>(n_pc);
      auto pc_coord_ = memManager_.template malloc<double>(3*n_pc);
      for(int i = 0; i < n_pc; i++){
        *(pc_+i) = double(aoints.molecule_.atoms[i].atomicNumber);
        for(int j = 0; j < 3; j++){
          *(pc_coord_+3*i+j) = aoints.molecule_.atoms[i].coord[j];
        }
      }
      pc = const_cast<double*>(pc_);
      pc_coord = const_cast<double*>(pc_coord_);
      result_ = efp_set_point_charges(efp_,n_pc,pc,pc_coord); 
     
    };

    // Compute one electron EFP contribution
    //
    // Coulomb part : monopole,dipole and quadrupole of EFP
    IntsT* One_electron_EFP_coulomb(void){
      auto efp_cou = aoints.computeEFPContributions_Coulomb(mult_coord, mult, &n_mult);
      return efp_cou;
    };
    
    // Polarization part
    IntsT* One_electron_EFP_pol(void){
      pol = memManager_.template malloc<double>(3*n_pol);
      pol_conj = memManager_.template malloc<double>(3*n_pol);
      efp_get_induced_dipole_values(efp_, pol);
      efp_get_induced_dipole_conj_values(efp_, pol_conj);
      for(int i = 0; i < n_pol*3; i++){
        *(pol+i) = ( ( *(pol+i) + *(pol_conj+i) ) / 2 );
      }
      auto efp_pol = aoints.computeEFPContributions_Pol(pol_coord, pol, &n_pol);
      memManager_.template free(pol);
      memManager_.template free(pol_conj);

      return efp_pol;
    };
 
    // Compute wavefunction dependent energy
    void Wavefunction_dependent(double* energy_ptr){
      result_ = efp_get_wavefunction_dependent_energy(efp_,energy_ptr);
    };



    // Begin EFP compute
    void EFP_Compute(int do_gradient){
      result_ = efp_compute(efp_, do_gradient);
    };

    ::efp_energy* EFP_get_contribution(){
      efp_get_energy(efp_,&efp_energy_);
      return &efp_energy_;
    };

    void get_new_user_data(SingleSlaterBase* ss_1){
      propagator_ = ss_1;
      efp_->get_electron_density_field_user_data = ss_1;
    }
      
 
    void Clean_up(void){
      efp_shutdown(efp_);
    };
   
    void Convert_to_string(void){
      efp_result_to_string(result_);
    };

  //
  //
  };


};
#endif
