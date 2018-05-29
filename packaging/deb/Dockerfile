FROM debian

WORKDIR /root

RUN apt --assume-yes update
RUN apt --assume-yes install python git
RUN apt --assume-yes install ruby ruby-dev rubygems build-essential
RUN gem install --no-ri --no-rdoc fpm
RUN git clone https://git.ligo.org/finesse/pykat.git
RUN apt --assume-yes install python-pip

CMD ["python", "/host/build.py"]