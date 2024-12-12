import sys
import openai

# Set OpenAI API key
openai.api_key = sys.argv[2]

# Get abstracts from command-line arguments
abstracts = sys.argv[1:]

# Function to use GPT-3 to extract disease-gene associations from abstracts
def extract_disease_gene_associations(abstracts):
    results = []
    for abstract in abstracts:
        response = openai.Completion.create(
            engine="text-davinci-003",  # Use GPT-3 engine
            prompt=f"Extract the disease-gene associations from the following scientific abstract: {abstract}",
            max_tokens=150,
            temperature=0.7
        )
        result = response.choices[0].text.strip()
        results.append(result)
    return results

# Extract disease-gene associations
disease_gene_associations = extract_disease_gene_associations(abstracts)

# Output results
print("\n".join(disease_gene_associations))
