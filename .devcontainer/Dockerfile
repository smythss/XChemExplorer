FROM gcr.io/diamond-privreg/xchem/ccp4:7.1 as ccp4
FROM gcr.io/diamond-privreg/xchem/phenix:1.20 as phenix

FROM rockylinux:9 as devel

COPY --from=ccp4 /ccp4-7.1 /ccp4-7.1
COPY --from=phenix /phenix-1.20 /phenix-1.20

RUN dnf update -y \
    && dnf install -y \
    libXrender fontconfig libXext \
    glib2 libSM libXi libXrandr libXfixes libXcursor libXinerama \
    libgomp libXdamage libxcb \
    python3-pip curl \
    && dnf clean all \
    && echo . /ccp4-7.1/bin/ccp4.setup-sh | tee ~/.bashrc ~/.zshrc \
    && curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o /tmp/get-pip.py \
    && /ccp4-7.1/bin/ccp4-python /tmp/get-pip.py \
    && /ccp4-7.1/bin/ccp4-python -m pip install flake8 \
    && pip3 install black[python2]

ENV QT_X11_NO_MITSHM=1 \
    PATH=$PATH:/ccp4-7.1/bin:/phenix-1.20/build/bin
