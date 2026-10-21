# Artifact image for "New IND-CPA-D Attacks on BFV and BGV".
#
# The three attacks need three different upstream builds, so each is built in
# its own stage and only the installed prefix is carried into the final image:
#
#   HElib   @ 3e337a6            Section 3, BGV bootstrapping
#   OpenFHE @ 2bf4524 (v1.3)     Section 5, BFV ct-ct, before the revised bound
#   OpenFHE @ df495ba (v1.5.0)   Section 4, BFV pt-ct, and Appendix E.3
#
# Build:  docker build -t indcpad-artifact .
# Run:    docker run -it --rm -v "$PWD/results:/artifact/results" indcpad-artifact

ARG BASE=ubuntu:24.04
# Compiling OpenFHE takes a few GB per job; raise only if the builder has the RAM.
ARG JOBS=4

############################  common build environment  ########################
FROM ${BASE} AS builder-base
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake git ca-certificates patch m4 \
      libgmp-dev libntl-dev libomp-dev \
 && rm -rf /var/lib/apt/lists/*

############################  HElib (BGV attack)  ##############################
FROM builder-base AS helib
ARG JOBS
ARG HELIB_COMMIT=3e337a6
RUN git clone https://github.com/homenc/HElib.git /src/HElib \
 && git -C /src/HElib checkout ${HELIB_COMMIT}
# The attack's coefficient extraction carries OpenMP pragmas that were being
# compiled away for want of -fopenmp. ENABLE_THREADS is a separate matter: it is
# what makes the collection's nthreads argument mean anything;
# without it HElib ignores the flag and runs on one core, which turns the day the
# paper reports into weeks. Ubuntu's NTL is built with NTL_THREADS, as it needs.
# The attack replaces a handful of HElib sources; they are overlaid before building.
COPY BGV/HElib-msater/include/helib/ /src/HElib/include/helib/
COPY BGV/HElib-msater/src/           /src/HElib/src/
COPY BGV/HElib-msater/tests/         /src/HElib/tests/
RUN cmake -S /src/HElib -B /src/HElib/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DGMP_DIR=/usr -DNTL_DIR=/usr \
      -DENABLE_TEST=ON -DPEDANTIC_BUILD=OFF \
      -DENABLE_THREADS=ON -DCMAKE_CXX_FLAGS=-fopenmp \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
 && cmake --build /src/HElib/build --target indcpad -j"${JOBS}" \
 && mkdir -p /opt/helib/bin && cp /src/HElib/build/bin/indcpad /opt/helib/bin/

# Both OpenFHE versions install a libOPENFHEpke.so.1, so LD_LIBRARY_PATH cannot
# serve them both; the final image gives each binary an RPATH to its own
# version's lib. It has to be done afterwards because OpenFHE's CMakeLists
# hardcodes set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE), overriding any -D flag.
############################  OpenFHE v1.3 (ct-ct, Section 5)  #################
FROM builder-base AS openfhe13
ARG JOBS
ARG OPENFHE_V13=2bf4524
RUN git clone https://github.com/openfheorg/openfhe-development.git /src/ofhe13 \
 && git -C /src/ofhe13 checkout ${OPENFHE_V13}
COPY BFV-ct-ct/bfv_ctct.hpp BFV-ct-ct/indcpad_ctct.cpp BFV-ct-ct/ctct_multikey.cpp \
     /src/ofhe13/src/pke/examples/
COPY ["BFV-ct-ct/(pke) CMakeLists.txt", "/src/ofhe13/src/pke/CMakeLists.txt"]
# The google-benchmark bundled with this commit does not configure with current
# toolchains, hence the two OFF switches. Dropping -Werror is for GCC 13, which
# is stricter about signed/unsigned comparison than the compiler the paper used;
# the attack sources are shipped exactly as they produced the published numbers.
RUN sed -i 's/-Wall -Werror/-Wall/' /src/ofhe13/CMakeLists.txt \
 && cmake -S /src/ofhe13 -B /src/ofhe13/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_BENCHMARKS=OFF -DBUILD_UNITTESTS=OFF \
 && cmake --build /src/ofhe13/build --target indcpad_ctct ctct_multikey -j"${JOBS}" \
 && mkdir -p /opt/openfhe-v13/bin /opt/openfhe-v13/lib \
 && cp -a /src/ofhe13/build/lib/libOPENFHE*.so* /opt/openfhe-v13/lib/ \
 && cp /src/ofhe13/build/bin/examples/pke/indcpad_ctct \
       /src/ofhe13/build/bin/examples/pke/ctct_multikey /opt/openfhe-v13/bin/

############################  OpenFHE v1.5.0 (pt-ct, Section 4, and Appendix E.3)
FROM builder-base AS openfhe15
ARG JOBS
ARG OPENFHE_V15=df495ba
RUN git clone https://github.com/openfheorg/openfhe-development.git /src/ofhe15 \
 && git -C /src/ofhe15 checkout ${OPENFHE_V15}
COPY BFV-pt-ct/bfv_ptct.hpp BFV-pt-ct/bfv_ptct_sk_general.cpp BFV-pt-ct/bfv_ptct_sk.cpp \
     BFV-pt-ct/bfv_ptct_param_prescreener.cpp \
     /src/ofhe15/src/pke/examples/
COPY BFV-ct-ct/bfv_ctct.hpp BFV-ct-ct/indcpad_ctct.cpp BFV-ct-ct/ctct_multikey.cpp \
     /src/ofhe15/src/pke/examples/
COPY ["BFV-pt-ct/(pke) CMakeLists.txt", "/src/ofhe15/src/pke/CMakeLists.txt"]
# The overlay carries a WITH_FIXED_SEED option, left OFF here. A fixed seed only
# makes runs identical within one build: OpenFHE draws the ternary secret through
# std::uniform_int_distribution, whose sequence the standard does not pin down, so
# libc++ and libstdc++ diverge from the same seed regardless. Every run therefore
# uses a fresh key, which is the more honest demonstration for an attack anyway.
COPY ["BFV-pt-ct/(main) CMakeLists.txt", "/src/ofhe15/CMakeLists.txt"]
COPY BFV-pt-ct/blake2engine.cpp /src/ofhe15/src/core/lib/utils/prng/blake2engine.cpp
RUN sed -i 's/-Wall -Werror/-Wall/' /src/ofhe15/CMakeLists.txt \
 && cmake -S /src/ofhe15 -B /src/ofhe15/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_BENCHMARKS=OFF -DBUILD_UNITTESTS=OFF \
 && cmake --build /src/ofhe15/build --target \
      bfv_ptct_sk_general indcpad_ctct ctct_multikey -j"${JOBS}" \
 && mkdir -p /opt/openfhe-v15/bin /opt/openfhe-v15/lib \
 && cp -a /src/ofhe15/build/lib/libOPENFHE*.so* /opt/openfhe-v15/lib/ \
 && cp /src/ofhe15/build/bin/examples/pke/bfv_ptct_sk_general \
       /src/ofhe15/build/bin/examples/pke/indcpad_ctct \
       /src/ofhe15/build/bin/examples/pke/ctct_multikey /opt/openfhe-v15/bin/

############################  final image  #####################################
FROM ${BASE}
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      libgmp10 libntl44 libgomp1 make patchelf \
      python3 python3-venv python3-pip \
 && rm -rf /var/lib/apt/lists/*

# The BGV analysis calls np.trapezoid, which is numpy 2.0 onward, whereas Ubuntu
# 24.04 ships 1.26; hence a virtual environment rather than the distribution
# packages. The versions are pinned so that a rebuild gets the same plots.
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir \
      numpy==2.2.6 scipy==1.15.3 pandas==2.2.3 matplotlib==3.10.3 seaborn==0.13.2
ENV PATH="/opt/venv/bin:${PATH}"

COPY --from=helib     /opt/helib        /opt/helib
COPY --from=openfhe13 /opt/openfhe-v13  /opt/openfhe-v13
COPY --from=openfhe15 /opt/openfhe-v15  /opt/openfhe-v15

# The libraries reference each other by the build-tree path too, so they need
# the same treatment; *.so and *.so.1 are symlinks and are left alone.
RUN for v in 13 15; do \
      patchelf --set-rpath /opt/openfhe-v$v/lib \
        /opt/openfhe-v$v/bin/* /opt/openfhe-v$v/lib/*.so.?.*; \
    done \
 && ! ldd /opt/openfhe-v13/bin/* /opt/openfhe-v15/bin/* | grep "not found"

WORKDIR /artifact
COPY artifact/Makefile artifact/README.md artifact/report.py artifact/to_npy.py ./
COPY BGV/BGV-BTS-KRD-slices-and-plot.py  BGV/
COPY BFV-ct-ct/ctct-boxplot.py           BFV-ct-ct/
COPY BFV-pt-ct/                          BFV-pt-ct/
COPY BFV-ct-ct/paper_results/           BFV-ct-ct/paper_results/

# There is no display here, so the scripts' plt.show() would drop the figures on
# the floor; the .pth below turns it into a savefig into FIGDIR.
COPY artifact/indcpad_figsave.py /opt/venv/lib/python3.12/site-packages/
RUN echo 'import indcpad_figsave' > /opt/venv/lib/python3.12/site-packages/zz-indcpad.pth
# Unbuffered, so that a long analysis shows progress rather than looking hung.
ENV MPLBACKEND=Agg FIGDIR=/artifact/results/figures PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/matplotlib

# The collected BGV data is large and ships beside the image; mount or copy it in.
VOLUME ["/artifact/data-collected-by-authors", "/artifact/results"]

ENV PATH="/opt/helib/bin:/opt/openfhe-v13/bin:/opt/openfhe-v15/bin:${PATH}"
CMD ["make", "help"]
