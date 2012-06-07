project(snputil)

find_package(Boost REQUIRED)

include_directories(${BOOST_INCLUDES})

add_executable(snpcord snpcord.cpp)
target_link_libraries(snpcord boost_program_options)

