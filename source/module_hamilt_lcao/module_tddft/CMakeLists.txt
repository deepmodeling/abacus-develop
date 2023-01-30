if(ENABLE_LCAO)
    list(APPEND objects
        ELEC_evolve.cpp
        LCAO_evolve.cpp
    )

    add_library(
        tddft
        OBJECT
        ${objects}
    )

    if(ENABLE_COVERAGE)
    add_coverage(tddft)
    endif()

endif()
