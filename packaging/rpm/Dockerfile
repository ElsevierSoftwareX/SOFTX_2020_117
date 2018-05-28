FROM scientificlinux/sl:7

WORKDIR /root

RUN yum -y update
RUN yum -y install git rpm-build gcc gsl-devel make gcc-c++
RUN git clone https://git.ligo.org/finesse/pykat.git

#CMD ["python", "/host/build.py"]