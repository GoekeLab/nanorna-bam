FROM nfcore/base
LABEL authors="Laura Wratten" \
      description="Docker image containing all requirements for nf-core/nanornabam pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml \
    && conda clean -a \
    && R -e "devtools::install_github('GoekeLab/bambu', dep = TRUE)" 
ENV PATH /opt/conda/envs/nf-core-nanornabam-1.0dev/bin:$PATH
