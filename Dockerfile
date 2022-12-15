# Python support can be specified down to the minor or micro version
# (e.g. 3.6 or 3.6.3).
# OS Support also exists for jessie & stretch (slim and full).
# See https://hub.docker.com/r/library/python/ for all supported Python
# tags from Docker Hub.
FROM continuumio/miniconda3:4.12.0

LABEL Name=tinycov Version=0.4.0

# Install python dependencies
COPY * ./ /app/
WORKDIR /app

RUN conda update -y conda
RUN conda config --add channels bioconda

# Get 3rd party packages directly from conda
RUN conda install -c conda-forge -y \
    samtools \
    htslib \
    pysam

RUN pip install -Ur requirements.txt
# Using pip:
RUN pip install .
#CMD ["python3", "-m", "hicstuff.main"]
ENTRYPOINT [ "tinycov" ]
