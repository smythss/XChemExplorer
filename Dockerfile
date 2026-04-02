FROM gcr.io/diamond-privreg/xchem/ccp4:7.1 as ccp4
FROM gcr.io/diamond-privreg/xchem/phenix:1.20 as phenix

FROM rockylinux:9 as xce

COPY --from=ccp4 /ccp4-7.1 /ccp4-7.1

ARG XCE_DIR=/xce
WORKDIR ${XCE_DIR}

COPY . ${XCE_DIR}

RUN dnf update -y \
    && dnf install -y \
        libXrender fontconfig libXext \
        glib2 libSM libXi libXrandr libXfixes libXcursor libXinerama \
        libgomp libXdamage libxcb \
    && dnf clean all

ENV QT_X11_NO_MITSHM=1 \
    XChemExplorer_DIR=${XCE_DIR}

ENTRYPOINT /ccp4-7.1/bin/ccp4-python -m xce
