# Use the official Python image as the base image
FROM python:3.9

# Set the working directory inside the container
WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the entire content of your Python project into the container's working directory
COPY . .

# Command to run your Python service
CMD ["python", "service_main.py"]

