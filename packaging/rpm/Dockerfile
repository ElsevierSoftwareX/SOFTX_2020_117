FROM scientificlinux/sl:7

WORKDIR /root

RUN yum -y update
RUN yum -y install git rpm-build ruby-devel gcc make rpm-build rubygems python python-setuptools
RUN git clone https://git.ligo.org/finesse/pykat.git
RUN gem install --no-ri --no-rdoc fpm

CMD ["python", "/host/build.py"]