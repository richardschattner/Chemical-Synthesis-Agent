# Chemical Synthesis Agent

An intelligent expert system for planning chemical synthesis pathways. This tool resolves natural language queries into chemical structures, searches a knowledge graph for valid reaction pathways, and generates step-by-step procedural guides using an LLM.

## 🧪 Overview

The **Chemical Synthesis Agent** is designed to solve the retrosynthesis problem—finding a pathway to synthesize a target molecule from a set of available starting materials. Unlike pure LLM approaches that may hallucinate non-existent reactions, this system grounds its reasoning in a verifiable **Neo4j Knowledge Graph** constructed from the **Open Reaction Database (ORD)**.

### How it Works
1.  **User Query**: The user asks, e.g., *"How can I synthesize Aspirin starting from Salicylic Acid and Acetic Anhydride?"*
2.  **Entity Resolution**: A **Parser Agent** extracts chemical entities and converts them to canonical SMILES strings using:
    *   **OPSIN API** (for systematic IUPAC names)
    *   **PubChem API** (for common/trade names)
    *   RDKit (for canonicalization)
3.  **Graph Search**: A Breadth-First Search (BFS) algorithm traverses the Neo4j database to find the shortest reaction path connecting the starting materials to the target product.
4.  **Data Retrieval**: Reaction conditions, procedures, catalysts, and solvents are fetched for every step in the found path.
5.  **Narration**: A **Narrator Agent** (LangChain + GPT-4o) synthesizes this structured data into a coherent, human-readable laboratory guide.

## 🚀 Features

*   **Natural Language to Graph**: Seamlessly bridges text queries with structured graph database queries.
*   **Exact Chemical Resolution**: Robust handling of chemical nomenclature (common names & IUPAC).
*   **Fact-Based Output**: All reactions are derived from actual experimental data (ORD), ensuring 0% hallucination on the *existence* of a reaction.
*   **Interactive Visualization**: A Streamlit frontend visualizes the reaction pathway using **PyVis**.
*   **Scalable Ingestion**: Optimized batch ingestion pipeline to load massive datasets (GBs of JSONL) into Neo4j.

## 🛠️ Tech Stack

*   **Orchestration**: [LangChain](https://www.langchain.com/) (LCEL)
*   **Database**: [Neo4j](https://neo4j.com/) (Graph Database)
*   **LLM**: OpenAI GPT-4o-mini
*   **Cheminformatics**: [RDKit](https://www.rdkit.org/)
*   **Frontend**: [Streamlit](https://streamlit.io/)
*   **APIs**: OPSIN, PubChem (PUG REST)

## 📦 Installation

### Prerequisites
*   Python 3.10+
*   Neo4j Database (Local Desktop or AuraDB)
*   OpenAI API Key

### Setup

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/richardschattner/Chemical-Synthesis-Agent.git
    cd Chemical-Synthesis-Agent
    ```

2.  **Install Dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

3.  **Configuration**:
    Create a `.env` file in the root directory:
    ```env
    OPENAI_API_KEY=sk-proj-...
    NEO4J_URI=bolt://localhost:7687
    NEO4J_USERNAME=neo4j
    NEO4J_PASSWORD=your_password
    ```

### Data Ingestion

1.  Download the **Open Reaction Database (ORD)** data files (JSONL format) and place them in the `data/` directory.
2.  Run the ingestion script to populate your Neo4j database:
    ```bash
    python ingest.py
    ```
    *Note: The script uses batch processing (`UNWIND`) to efficiently handle millions of nodes.*

## 🖥️ Usage

**Run the Streamlit Application**:
```bash
streamlit run app.py
```

### Example Queries
*   *"Synthesize paracetamol from 4-aminophenol and acetic anhydride."*
*   *"Find a path to make Benzocaine using p-aminobenzoic acid and ethanol."*
*   *"How do I make Aspirin?"* (The system will infer common starting materials if not specified, though explicit inputs improve precision).

## 📂 Project Structure

*   `app.py`: Main Streamlit application entry point.
*   `agent.py`: LangChain orchestration logic.
*   `retro.py`: Core BFS retrosynthesis algorithms and Graph DB queries.
*   `ingest.py`: ETL pipeline for loading ORD data into Neo4j.
*   `parser_agent.py`: Agent responsible for extracting and resolving chemical entities.
*   `narrator_agent.py`: Agent responsible for generating natural language explanations.

## 🛡️ License

This project is open-source and available under the [MIT License](LICENSE).

---
*Note: This system relies on the quality and coverage of the data loaded into Neo4j. If a reaction path exists in chemical literature but is not in the ORD subset loaded, the agent will not find it.*
