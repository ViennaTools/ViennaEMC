cmake_minimum_required(VERSION 3.4)

project("mosfet2D")

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNAEMC_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${VIENNAEMC_LIBRARIES})

# Add to example target
add_dependencies(buildExamples ${PROJECT_NAME})