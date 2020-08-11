# Copyright (c) 2020 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

FROM diehlpk/nonlocalmodels:baseimage

WORKDIR /app

RUN git clone  https://github.com/nonlocalmodels/nonlocalheatequation.git

RUN cd /app/nonlocalheatequation && mkdir build && cd build && cmake .. -DHPX_IGNORE_COMPILER_COMPATIBILITY=ON -DCMAKE_BUILD_TYPE=Release && make -j && sudo make install

