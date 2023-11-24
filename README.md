# Gap-filling workflow for metabolic networks with BioMeneco

This project is an integrated workflow combining the capabilities of tools like BioISO, Meneco and other relevant developed methods to address gaps in Genome-scale Metabolic models (GEMs). These gaps often arise from database limitations and incorrect genome annotations. While existing tools that identify network deficiencies or suggest reactions to fill these gaps, their efficiency decreases with increasing network complexity. This workflow aims to streamline this process, offering a more efficient and automated approach

## Installation
### Manual

#### Create a conda environment with the following command

```bash
conda env create -n gapfilling
```

#### Activate the environment with the following command

```bash
conda activate gapfilling
```

#### Clone the repository
```bash
git clone https://github.com/BioSystemsUM/biomeneco/tree/gap_filler_model
```

#### Navigate to the repository folder
```bash
cd [repository folder name]
# Replace [repository folder name] with the name of the folder created by cloning the repository
```

#### Install requirements with the following command

```bash
pip install -r requirements.txt
```


### Docker Installation and Setup

#### Install Docker
First, ensure Docker is installed on your system. If it's not installed, you can download it from the [Docker website](https://www.docker.com/).

#### Clone the Repository
Similar to the manual installation, clone the GitHub repository to your local machine.

```bash
git clone https://github.com/BioSystemsUM/biomeneco.git
```

#### Navigate to the Repository Folder
Navigate to the folder where the Dockerfile is present. Replace [biomeneco] with the actual name of the folder created by cloning the repository.

```bash
cd biomeneco
```

### Build the Docker Image
Build the Docker image from the Dockerfile in the repository. This step assumes you are in the root directory of the cloned repository.

```bash
   docker build . -t gapfilling
```



## Usage

### Manual Usage

#### Data Preparation
Place your `model.xml` (rename your model to be gap-filled to 'model.xml') and a JSON file with submission parameters in the input folder. The JSON file should follow the structure of the example provided in the repository, available here:
[ExampleSubmissionParameters.JSON](https://github.com/BioSystemsUM/biomeneco/blob/gap_filler_model/ExampleSubmissionParameters.JSON)

#### Explaining User Inputs
- **SBML Model:** The draft metabolic model submitted by the user must be in SBML format and renamed to `model.xml`. This model is read by [COBRApy](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-7-74) to allow easy manipulation.
- **JSON Parameters File:** This file contains key parameters for processing the metabolic model. These include:
  1. `objective_function_id` (mandatory): Specifies the objective function of the metabolic model for gap-filling (e.g., `e_Biomass__cytop`).
  2. `sinks` (optional): Lists metabolites used as seeds. For each metabolite, a sink reaction is generated.
  3. `compartments` (optional): Enumerates different cellular compartments in the model, defining which compartments will be replicated in the universal model.
  4. `max_solutions` (optional): Sets the maximum number of solution sets the algorithm should consider.

#### Running the Workflow
Execute the workflow using `MenecoRunner.py`. This script requires two arguments:
- `input_folder`: The path to the input folder containing `model.xml` and the JSON file described above.
- `output_folder`: The path to the output folder where results will be stored.

```bash
python MenecoRunner.py /path/to/input_folder /path/to/output_folder
```


### Docker Usage

Using Docker simplifies the setup process. Follow these steps to run the workflow in a Docker container:

1. **Prepare Your Data:**
   - Rename your SBML model file to `model.xml` and place it in a designated input folder on your local machine. 
   - Include a JSON file with submission parameters (`SubmissionParameters.json`) in the same input folder. The structure of this file should follow the example provided in the repository, which can be found here: [ExampleSubmissionParameters.JSON](https://github.com/BioSystemsUM/biomeneco/blob/gap_filler_model/ExampleSubmissionParameters.JSON).

2. **Run the Docker Container:**
   - Start the Docker container, specifying the paths to your input and output folders. The container will automatically execute `MenecoRunner.py` using the provided directories.
   - Replace `your_path_to_package/gap_filling_dl/`, `/path/to/input_folder`, and `/path/to/output_folder` with the actual paths to your input and output directories.

   ```bash
   docker run -v your_path_to_package/gap_filling_dl/:/home gapfilling python3 /home/src/gap_filling_dl/MenecoRunner.py /path/to/input_folder /path/to/output_folder
    ```


## Understanding the Results and Report Structure

After running the gap-filling workflow, a comprehensive report is generated in JSON format. This report provides a detailed overview of the process, including the files used, execution time, and the specific reactions and metabolites involved. Below is a breakdown of the report structure and an example of the final results.

### Report Structure

The report includes the following sections:

- **Execution Time**: The total time taken for the algorithm to execute, measured in seconds.
- **Objective**: The objective function value achieved by the model, after optimization with COBRApy, referring to the growth rate.
- **Artificial Removed Reactions**: Reactions that were artificially removed from the model, with their identifiers. (Note: This is for test usage only and may not be relevant in all cases.)
- **Unproducible Targets**: Targets that could not be produced by the model, listed with their identifiers.
- **Unreconstructable Targets**: Targets that could not be reconstructed by the model.
- **Reconstructable Targets**: Targets that were successfully reconstructed, along with their identifiers.
- **Minimal Completion**: A list of reactions considered for the minimal completion of the network.
- **Additional Seeds**: Supplementary seeds identified during the process, after their initial recognition using the `Model` class.
- **Essential Additional Seeds**: Seeds deemed crucial for the metabolic network’s functionality, identified after the initial seeds have been added.
- **Additional Demands**: Lists the demand reactions added to the model to address dead ends.
- **Summary**: Provides an overview of the filled model, detailing reactions, associated metabolites, factors, and flux values.
- **Positive Solutions**: Solutions that successfully produce biomass, indicating viable modifications to the metabolic model.












