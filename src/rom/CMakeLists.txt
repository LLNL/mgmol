set(SOURCES 
  romOffline.cc 
)

add_library(mgmol_rom ${SOURCES})
include_directories("${LIBROM_PATH}/lib")
install(TARGETS mgmol_rom DESTINATION lib)
