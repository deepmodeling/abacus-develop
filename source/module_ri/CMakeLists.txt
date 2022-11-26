if (ENABLE_LIBRI)
    add_library(
            libri
            OBJECT
            Matrix_Orbs11.cpp
            Matrix_Orbs21.cpp
            RI_2D_Comm.cpp
    )
endif()