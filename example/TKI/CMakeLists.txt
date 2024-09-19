include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(tki_general STATIC tki_general.cxx)


add_library(nuwroTKI SHARED nuwroTKI.cxx)
add_library(genieTKI SHARED genieTKI.cxx)

target_link_libraries(nuwroTKI tki_general ROOT::Physics ROOT::ROOTDataFrame)
target_link_libraries(genieTKI tki_general genie_general ROOT::Physics ROOT::ROOTDataFrame)
