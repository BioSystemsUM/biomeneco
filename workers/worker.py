from flask import Flask, request, send_file, jsonify
from logging.handlers import TimedRotatingFileHandler
import subprocess
import os
import shutil
import datetime
import logging
import time
import re

from FilesUtilities import compressFiles

app = Flask(__name__)

logPath = '/workdir/logs/'

if not os.path.exists(logPath):
    os.makedirs(logPath)

# format the log entries
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
handler = TimedRotatingFileHandler(logPath + 'MenecoWorker.log', when='midnight', backupCount=20)
handler.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

PROCESS_PATH_MENECO = '/home/MenecoRunner.py'
SUBMISSIONS_PATH = '/workdir/workerSubmissions/'
RESULTS_PATH = '/workdir/resultsWorker/'
OBJECTIVE = 'e_Biomass__cytop'

if not os.path.exists(SUBMISSIONS_PATH):
    os.makedirs(SUBMISSIONS_PATH)
if not os.path.exists(RESULTS_PATH):
    os.makedirs(RESULTS_PATH)

start_time = datetime.datetime.now()

running = False

logger.info("Meneco Worker running...")


@app.route("/Meneco", methods=["POST"])
def startMeneco():
    global start_time
    global running

    logger.info("Submission in progress...")

    start_time = datetime.datetime.now()

    if os.path.exists(SUBMISSIONS_PATH):
        shutil.rmtree(SUBMISSIONS_PATH, ignore_errors=True)

    if os.path.exists(RESULTS_PATH):
        shutil.rmtree(RESULTS_PATH, ignore_errors=True)

    os.makedirs(SUBMISSIONS_PATH)
    os.makedirs(RESULTS_PATH)

    logger.debug("Directories reset complete")

    try:
        for file in request.files.values():
            destination = SUBMISSIONS_PATH + str(file.filename)
            logger.info("Saving file: " + destination)
            file.save(destination)

        subprocess.Popen(["python3", PROCESS_PATH_MENECO, SUBMISSIONS_PATH, RESULTS_PATH, OBJECTIVE])

        running = True

        logger.info("returning code 102!")

        return ('processing', 102)

    except IndexError:
        logger.info('An error occurred while processing the submission, returning code 500!')
        return ('error', 500)


@app.route('/status')
def display_msg():
    global running

    files = os.listdir(RESULTS_PATH)

    if not running:
        return jsonify({"result": "Not running"}), 410

    # if abs(datetime.datetime.now() - start_time) > datetime.timedelta(hours=2):
    #    return jsonify({"result": "Time out"}), 408

    if "processComplete" in files:

        if "400" in files:
            with open(RESULTS_PATH + "/400", "r") as error_file:
                message = str(error_file.readlines()[0])
            return jsonify({"result": message}), 400

        return send_file(RESULTS_PATH + "results.zip", as_attachment=True,
                         download_name='results.zip'), 200
    #
    # StartTime = 0
    # TotalSequences = 0
    # ProcessedSequences = 0
    #
    # if "StartTime" in os.listdir(SUBMISSIONS_PATH):
    #     with open(SUBMISSIONS_PATH + "StartTime", "r") as StartTimeFile:
    #         StartTime = float(StartTimeFile.readline())
    #
    # if "TotalSequences" in os.listdir(SUBMISSIONS_PATH):
    #     with open(SUBMISSIONS_PATH + "TotalSequences", "r") as TotalSequencesFile:
    #         TotalSequences = int(TotalSequencesFile.readline())
    #
    # if "ProcessedSequences" in os.listdir(SUBMISSIONS_PATH):
    #     with open(SUBMISSIONS_PATH + "ProcessedSequences", "r") as ProcessedSequencesFile:
    #         ProcessedSequences = int(ProcessedSequencesFile.readline())
    #
    # AverageTime = float(time.monotonic()-StartTime)
    #
    # if ProcessedSequences > 0:
    #     AverageTime = float((time.monotonic()-StartTime)/ProcessedSequences)
    #
    # ExpectedTime = float((TotalSequences - ProcessedSequences)*AverageTime)

    time = re.compile("\[.*\]")

    indexation_time = 0.0
    iteration_time = 0.0
    start = False
    iteration_end = False

    total_blocks = 1
    total_shapes = 1
    total_chunks = 1

    counter = 0

    processed_blocks = 0
    processed_shapes = 0
    processed_chunks = 0

    expected_time = 0.0

    if os.path.exists("/workdir/diamond.log"):

        shutil.copy("/workdir/diamond.log", "/workdir/diamond-temp.log")

        with open("/workdir/diamond-temp.log", "r") as logger:

            for line in logger.readlines():

                if not start:
                    for match in re.findall(time, line):
                        indexation_time += float(match[1:-2])

                    if line.startswith("Processing query block"):
                        line = line.replace("\n", "").replace(".", "")
                        for line_part in line.split(","):
                            if "reference block" in line_part:
                                line_part = str(line_part).replace("reference block", "").replace(" ", "")
                                processed_blocks, total_blocks = line_part.split("/")
                            elif "shape" in line_part:
                                line_part = str(line_part).replace("shape", "").replace(" ", "")
                                processed_shapes, total_shapes = line_part.split("/")
                            elif "index chunk" in line_part:
                                line_part = str(line_part).replace("index chunk", "").replace(" ", "")
                                processed_chunks, total_chunks = line_part.split("/")

                        start = True

                else:

                    if line.startswith("Processing query block"):
                        counter += 1

                    if not iteration_end:

                        for match in re.findall(time, line):
                            iteration_time += float(match[1:-2])

                        if line.startswith("Processing query block"):
                            iteration_end = True

        os.remove("/workdir/diamond-temp.log")

    if start:
        expected_time = float(iteration_time) * int((int(total_blocks) * int(total_shapes) * int(total_chunks)) - int(counter))
    else:
        expected_time = -1.0

    return jsonify({"ProcessedSequences": str(counter),
                    "TotalSequences": str(int(total_blocks) * int(total_shapes) * int(total_chunks)),
                    "ExpectedTime": str(expected_time),
                    "result": "running"}), 202


@app.route('/handshake')
def handshake():
    logger.info("Handshake requested!")

    return jsonify({"result": "alive"}), 200


@app.route('/runningStatus')
def runningStatus():
    if os.path.exists(RESULTS_PATH):

        files = os.listdir(RESULTS_PATH)

        if "processComplete" in files or len(files) == 0:
            return jsonify({"result": "not running"}), 200

    return jsonify({"result": "running"}), 202


@app.route('/retrieveLogs')
def retrieveLogs():
    logsPath = '/workdir/logs'
    logsZip = '/logs.zip'

    if os.path.exists(logsZip):
        os.remove(logsZip)

    compressFiles(logsPath, logsZip)

    if not os.path.exists(logsZip):
        return jsonify({"Result": "Request unavailable!"}), 503

    return send_file(logsZip, as_attachment=True, download_name='logs.zip')


@app.route('/checkIfAnySubmissionWasBrutallyInterrupted')
def checkIfAnySubmissionWasBrutallyInterrupted():
    files = os.listdir(RESULTS_PATH)
    if "processComplete" not in files:
        with open(RESULTS_PATH + "processComplete", "w") as processCompleteError:
            processCompleteError.write("There has been an issue running the last submission and it was brutally interrupted.")

        return jsonify({"result": "There has been an issue running the last submission and it was brutally interrupted."}), 200

    return jsonify({"result": "No submission was brutally interrupted."}), 200


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=80, threaded=False, debug=True)
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
