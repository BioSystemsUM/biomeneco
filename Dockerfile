# Use the official Python image as the base image
FROM ubuntu:22.04

COPY ./ /home

RUN apt-get update
RUN apt-get install -y python3-pip python3.9
RUN apt-get install wget -y
RUN apt-get install -y git

# Install Python dependencies
RUN pip3 install Flask
RUN pip3 install --no-cache-dir -r /home/requirements.txt
RUN git clone https://github.com/BioSystemsUM/BioISO.git && cd BioISO && python3 setup.py install

RUN sh /home/utilities/cplex_studio2210.linux_x86_64.bin -f /home/utilities/cplex_response.txt

ENV PYTHONPATH "${PYTHONPATH}:/home/src/"

# Set the working directory inside the container
WORKDIR /workdir



EXPOSE 80

# Command to run your Python service
#CMD ["python3", "/home/src/gap_filling_dl/MenecoRunner.py"]

