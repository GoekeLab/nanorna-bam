FROM nfcore/base
LABEL authors="Laura Wratten" \
      description="Docker image containing all requirements for nf-core/nanornabam pipeline"

COPY environment.yml /


RUN conda env create -f /environment.yml \
    && conda clean -a

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
	&& rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
        && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

ENV R_BASE_VERSION 4.0.1
    
RUN apt-get update \
	&& apt-get install -t unstable -y --no-install-recommends \
	              aptitude \
	              libcurl4-openssl-dev \
	              libxml2-dev \
                gcc-9-base \
                r-base=${R_BASE_VERSION}-* \
		            r-base-dev=${R_BASE_VERSION}-* \
		            r-recommended=${R_BASE_VERSION}-* \
		            r-cran-xml \
                r-cran-devtools \
                r-cran-rcurl\
  && R -e "devtools::install_github('Goekelab/bambu')" \
  && rm -rf /tmp/downloaded_packages/ 
     



ENV PATH /opt/conda/envs/nf-core-nanornabam-1.0dev/bin:$PATH
