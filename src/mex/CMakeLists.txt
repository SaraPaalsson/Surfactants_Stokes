include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)


## MEX functions
matlab_add_mex(
	NAME mex_applymodifs
	SRC mex_applymodifs.cpp
)


matlab_add_mex(
	NAME mex_dospecquad3
	SRC mex_dospecquad3.cpp
)

matlab_add_mex(
	NAME mex_gmresinner_el
	SRC mex_gmresinner_el.c
)

matlab_add_mex(
	NAME mex_M1
	SRC mex_M1.c
)

matlab_add_mex(
	NAME mex_M1M4
	SRC mex_M1M4.c
)

matlab_add_mex(
	NAME mex_M1M5
	SRC mex_M1M5.c
)

matlab_add_mex(
	NAME mex_Trap2GL
	SRC mex_Trap2GL.c
)

matlab_add_mex(
	NAME mex_Trap2GL_real
	SRC mex_Trap2GL_real.c
)










