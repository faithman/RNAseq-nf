FROM continuumio/miniconda:latest

RUN conda config --add channels bioconda && \
     conda config --add channels conda-forge && \
     conda config --add channels defaults
          
RUN conda install kallisto=0.43.1 \
                  fastqc=0.11.7 \
                  multiqc=1.5 \
                  salmon=0.9.1 \
                  r-sleuth=0.29.0 




