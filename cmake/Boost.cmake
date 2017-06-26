set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED system)

include_directories(${Boost_INCLUDE_DIRS})
list(APPEND CQ_EXT_LINK ${Boost_LIBRARIES})
