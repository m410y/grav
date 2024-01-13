cmake_minimum_required(VERSION 2.8...3.28)
project(grav)

file(GLOB_RECURSE sources src/main/*.c src/main/*.h)
file(GLOB_RECURSE sources_test src/test/*.c)
file(GLOB_RECURSE data resources/*)

add_executable(particle ${sources} ${data})
target_compile_options(particle PUBLIC -Wall -Wextra -g)
target_include_directories(particle PUBLIC src/include)

target_link_libraries(particle PUBLIC
    m
)

find_package(Check)

if (Check_FOUND)
    add_executable(unit_tests ${sources_test} ${sources})
    target_compile_definitions(unit_tests PUBLIC UNIT_TESTS)
    target_compile_options(unit_tests PUBLIC -Wall -Wextra -g)
    target_include_directories(unit_tests PUBLIC src/test)
    set_target_properties(particle PROPERTIES ENABLE_EXPORTS on)
    target_link_libraries(unit_tests PUBLIC
        check
        particle
    )
endif()
