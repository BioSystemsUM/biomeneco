FROM rcarmo/ubuntu-python

USER root

#update linux and accessories
RUN apt-get -y update
RUN apt-get -y upgrade

#get pip
RUN apt-get install -y python-dev build-essential
#RUN apt-get install python3-pip

COPY requirements.txt ./home
RUN pip3 install --upgrade pip
RUN pip3 install -r /home/requirements.txt
ENV PYTHONPATH "/home:${PYTHONPATH}"

CMD mkdir /home/sharing_point
CMD mkdir /menecopy

ADD tests/run_bioiso_and_meneco.py /home/
ADD src/gap_filling_dl/biomeneco/BioISO/menecopy /home/menecopy
ADD tests/data/models /home/models


WORKDIR home



#CMD ["python","cplex/python/setup.py","install"]
CMD ["python","run_bioiso_and_meneco.py"]