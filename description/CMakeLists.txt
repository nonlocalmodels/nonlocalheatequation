#  Copyright (c) 2020 Patrick Diehl
#
#  SPDX-License-Identifier: BSL-1.0
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

find_package(LATEX COMPONENTS PDFLATEX)

include(UseLATEX.cmake)

add_latex_document(problem_description.tex IMAGE_DIRS ./imgs  FORCE_PDF)

add_custom_command(TARGET pdf POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy problem_description.pdf ${PROJECT_SOURCE_DIR}/description
)
