
# **PubMed-GPT3-Literature-Mining**

## Full Workflow: Mining Literature from PubMed Using PubMedGPT (GPT-3)

This project implements a full workflow to mine scientific literature from PubMed, utilizing GPT-3 to extract disease-gene associations. The workflow involves querying PubMed for relevant articles based on a search term, fetching their abstracts, and using GPT-3 to extract associations between diseases and genes. The results are then saved to a CSV file.

## **Directory Structure:**

```
/PubMed-GPT3-Literature-Mining
├── /WDL
│   ├── pubmed_gpt3.wdl                # WDL workflow file
│   ├── inputs.json                    # input file
│   ├── extract_disease_gene_associations.py # Script to extract disease-gene associations using GPT-3
│   ├── fetch_abstracts.py            # Script to fetch abstracts from PubMed based on article IDs
│   ├── fetch_pubmed_articles.py     # Script to fetch PubMed article IDs based on a search term
│   ├── save_results_to_csv.py       # Script to save extracted results to a CSV file
├── pubmed_gpt3_mining.py            # Main script to run the entire PubMed-GPT3 mining process
├── requirements.txt                # List of required Python packages
└── README.md                       # Project documentation
```

## **Project Overview:**

This project is designed to automate the process of literature mining for disease-gene associations using data from PubMed. The workflow leverages the **BioPython** library to query PubMed, **OpenAI's GPT-3** for natural language processing, and **Cromwell** to run the workflow defined in **WDL** (Workflow Description Language).

The workflow consists of the following tasks:
1. **Fetching PubMed article IDs** based on a search term.
2. **Fetching abstracts** from PubMed for the fetched article IDs.
3. **Using GPT-3** to extract disease-gene associations from the abstracts.
4. **Saving the results** to a CSV file for further analysis.

## **Installation:**

To get started, clone the repository and set up the environment:

```bash
git clone https://github.com/your-username/PubMed-GPT3-Literature-Mining.git
cd PubMed-GPT3-Literature-Mining
```

### **Install Required Python Packages:**

Create and activate a Python virtual environment (optional but recommended), then install the required packages:

```bash
# Create a virtual environment (optional)
python3 -m venv venv
source venv/bin/activate   # On Windows, use `venv\Scripts\activate`

# Install required Python packages
pip install -r requirements.txt
```

### **Requirements:**
- Python 3.x
- `BioPython` for PubMed access
- `OpenAI` for GPT-3 API
- `requests` for API interactions

You can view the dependencies in `requirements.txt`:

```text
biopython
openai
requests
```

## **Python Scripts:**

### **1. `fetch_pubmed_articles.py`:**
Fetches PubMed article IDs based on a search term using the Entrez API.

**Usage:**
```bash
python3 fetch_pubmed_articles.py --search_term "coronary artery disease" --max_results 10 --email "your-email@example.com"
```

### **2. `fetch_abstracts.py`:**
Fetches abstracts from PubMed based on the article IDs.

**Usage:**
```bash
python3 fetch_abstracts.py --article_ids "12345678" "23456789" "34567890"
```

### **3. `extract_disease_gene_associations.py`:**
Uses GPT-3 to extract disease-gene associations from the given abstracts.

**Usage:**
```bash
python3 extract_disease_gene_associations.py --abstracts "Abstract 1" "Abstract 2" --openai_api_key "your-openai-api-key"
```

### **4. `save_results_to_csv.py`:**
Saves the extracted disease-gene associations into a CSV file.

**Usage:**
```bash
python3 save_results_to_csv.py --disease_gene_associations "Association 1" "Association 2"
```

### **5. `pubmed_gpt3_mining.py`:**
Main Python script to execute the entire PubMed-GPT3 mining process.

This script orchestrates the entire process by calling the individual scripts for each task and processing the results. You can specify the search term, number of results, and OpenAI API key here.

**Usage:**
```bash
python3 pubmed_gpt3_mining.py --search_term "cancer" --max_results 10 --openai_api_key "your-openai-api-key" --email "your-email@example.com"
```

---

## **WDL Workflow:**

### **Workflow: `pubmed_gpt3.wdl`**

The workflow in **WDL** coordinates all tasks:
- **`fetch_pubmed_articles`**: Fetches article IDs from PubMed.
- **`fetch_abstracts`**: Fetches abstracts from PubMed.
- **`extract_disease_gene_associations`**: Uses GPT-3 to extract associations.
- **`save_results_to_csv`**: Saves the results to a CSV file.

To run this WDL workflow, you can use **Cromwell** or another WDL execution engine.

**Example JSON input file: `input.json`**:
```json
{
  "search_term": "Alzheimer's Disease",
  "max_results": 10,
  "openai_api_key": "your-openai-api-key",
  "email": "your-email@example.com"
}
```

Run the workflow using **Cromwell**:

```bash
java -jar cromwell.jar run pubmed_gpt3.wdl --inputs input.json
```

### **Output:**
The output will be a CSV file containing the disease-gene associations, saved as `disease_gene_associations.csv`.

---

## **Contributing:**

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make your changes and commit them (`git commit -am 'Add new feature'`).
4. Push to the branch (`git push origin feature-branch`).
5. Open a pull request to the main repository.

---

## **License:**

This project is licensed under the MIT License.

---

This `README.md` provides detailed instructions on setting up and using the project. It includes descriptions of the main components, installation steps, usage instructions for each Python script, and how to run the WDL workflow using Cromwell.
