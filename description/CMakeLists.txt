find_package(LATEX COMPONENTS PDFLATEX)

include(UseLATEX.cmake)

add_latex_document(problem_description.tex IMAGE_DIRS ./imgs  FORCE_PDF)

add_custom_command(TARGET pdf POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy problem_description.pdf ${PROJECT_SOURCE_DIR}/description
)
