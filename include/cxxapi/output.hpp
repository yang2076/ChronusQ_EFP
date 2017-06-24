#ifndef _INCLUDED_OUTPUT_HPP_
#define _INCLUDED_OUTPUT_HPP_


#include <chronusq_sys.hpp>

namespace ChronusQ {

  // String constants

  constexpr char  bannerTop[] = "--------------------------------------------------------------------------------";
  constexpr char  bannerMid[] = "  ----------------------------------------------------------------------------  ";
  constexpr char  bannerEnd[] = "--------------------------------------------------------------------------------";
  constexpr char  BannerTop[] = "================================================================================";
  constexpr char  BannerMid[] = "  ============================================================================  ";
  constexpr char  BannerEnd[] = "================================================================================";

  constexpr char CQBanner[] = 
"    ______ __                                      ____  \n "
"   / ____// /_   _____ ____   ____   __  __ _____ / __ \\ \n "
"  / /    / __ \\ / ___// __ \\ / __ \\ / / / // ___// / / / \n "
" / /___ / / / // /   / /_/ // / / // /_/ /(__  )/ /_/ /  \n "
" \\____//_/ /_//_/    \\____//_/ /_/ \\__,_//____/ \\___\\_\\  \n ";


  /**
   *  Dump the contents of a file to a specified output device.
   *
   *  \param [in]  out    Output device
   *  \param [out] fName  Name of file to output.
   */ 
  inline void CatFile(std::ostream &out, std::string fName) {

    std::ifstream inFile(fName);
    std::string line;
    while(std::getline(inFile,line))
      out << line << std::endl;
    
  }; // CatFile

  /**
   *  Ouput the standard ChronusQ file header to a specified 
   *  output defice.
   *
   *  \param [in]  out    Output device
   */ 
  inline void CQOutputHeader(std::ostream &out) {

    time_t currentClockTime;
    time(&currentClockTime);

    out << "ChronusQ Job Started: " << ctime(&currentClockTime) << std::endl;
    out << CQBanner << std::endl;
    
    out << "Release Version: "  
       << ChronusQ_VERSION_MAJOR << "." 
       << ChronusQ_VERSION_MINOR << "."
       << ChronusQ_VERSION_PATCH << std::endl;

    out << std::endl << std::endl << "Contributors List:" << 
        std::endl << BannerTop << std::endl;
    CatFile(out,ChronusQ_AUTHOR_LIST);

  }; // CQOuputHeader

  /**
   *  Ouput the standard ChronusQ file footer to a specified 
   *  output defice.
   *
   *  \param [in]  out    Output device
   */ 
  inline void CQOutputFooter(std::ostream &out) {

    time_t currentClockTime;
    time(&currentClockTime);

    out << std::endl << std::endl;
    out << "ChronusQ Job Ended: " << ctime(&currentClockTime) << std::endl;

  }; // CQOutputFooter

};


#endif
