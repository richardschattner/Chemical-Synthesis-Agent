# Chemical Synthesis Agent

An intelligent expert system for planning chemical synthesis pathways. This tool resolves natural language queries into chemical structures, searches a knowledge graph for valid reaction pathways, and generates step-by-step procedural guides using an LLM.
For an explanation of the architecture, algorithm and agent. Watch the video below

<p align="center">
  <a href="https://www.youtube.com/watch?v=dGVC0PSJO2A">
    <img src="https://img.youtube.com/vi/dGVC0PSJO2A/maxresdefault.jpg" width="600">
  </a>
</p>


## Overview

The **Chemical Synthesis Agent** is designed to solve the retrosynthesis problem—finding a pathway to synthesize a target molecule from a set of available starting materials. Unlike pure LLM approaches that may hallucinate non-existent reactions, this system grounds its reasoning in a verifiable **Neo4j Knowledge Graph** constructed from the **Open Reaction Database (ORD)**.

### How it Works
1.  **User Query**: The user asks, e.g., *"How can I synthesize Aspirin starting from Salicylic Acid and Acetic Anhydride?"*
2.  **Entity Resolution**: A **Parser Agent** extracts chemical entities and converts them to canonical SMILES strings using:
    *   **OPSIN API** (for systematic IUPAC names)
    *   **PubChem API** (for common/trade names)
    *   RDKit (for canonicalization)
3.  **Graph Search**: A Breadth-First Search (BFS) algorithm traverses the Neo4j knowledge graph to find the shortest reaction path connecting the starting materials to the target product.
4.  **Data Retrieval**: Reaction conditions, procedures, catalysts, and solvents are fetched for every step in the found path.
5.  **Narration**: A **Narrator Agent** (LangChain + GPT-4o) synthesizes this structured data into a coherent, human-readable laboratory guide.

## Features

*   **Natural Language to Graph**: Seamlessly bridges text queries with structured graph database queries.
*   **Exact Chemical Resolution**: Robust handling of chemical nomenclature (common names & IUPAC).
*   **Fact-Based Output**: All reactions are derived from actual experimental data (ORD), ensuring no hallucination on the *existence* of a reaction.
*   **Interactive Visualization**: A Streamlit frontend visualizes the reaction pathway using **PyVis**.
*   **Scalable Ingestion**: Optimized batch ingestion pipeline to load massive datasets (GBs of JSONL) into Neo4j.

## Tech Stack

*   **Orchestration**: [LangChain]
*   **Database**: [Neo4j] (Graph Database)
*   **LLM**: OpenAI GPT-4o-mini
*   **Cheminformatics**: [RDKit]
*   **Frontend**: [Streamlit]
*   **APIs**: OPSIN, PubChem (For converting names of chemicals into SMILES codes)

## Installation

### Prerequisites
*   Python 3.10+
*   Neo4j Database (Local Desktop or AuraDB)
*   OpenAI API Key (or any other LLM API provider)

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

1.  Download the **Open Reaction Database (ORD)** data files (pb.gz format)
2.  Convert the files into JSONL Format using the conversion script
    ```bash
    python convert.py
    ```
3.  Enforce uniqueness constraints for Molecules and Reactions in Neo4j to avoid duplicates by pasting the following into the Neo4j Console:
    ```bash
    CREATE CONSTRAINT molecule_smiles IF NOT EXISTS
    FOR (m:Molecule) REQUIRE m.smiles IS UNIQUE;
   
    CREATE CONSTRAINT reaction_id IF NOT EXISTS
    FOR (r:Reaction) REQUIRE r.reaction_id IS UNIQUE;
    ```
4.  Run the ingestion script to populate your Neo4j database:
    ```bash
    python ingest.py
    ```

## Usage

**Run the Streamlit Application**:
```bash
streamlit run app.py
```

### Example Queries
*   *"Synthesize paracetamol from 4-aminophenol and acetic anhydride."*
*   *"Find a path to make Benzocaine using p-aminobenzoic acid and ethanol."*
*   *"How can I synthesize aspirin from salicylic acid and acetic anhydride?"*
  
![The results of an example query](/example_output.png)

## Project Structure

*   `app.py`: Main Streamlit application entry point.
*   `agent.py`: LangChain orchestration logic.
*   `retro.py`: Core BFS retrosynthesis algorithms and Graph DB queries.
*   `ingest.py`: ETL pipeline for loading ORD data into Neo4j.
*   `parser_agent.py`: Agent responsible for extracting and resolving chemical entities.
*   `narrator_agent.py`: Agent responsible for generating natural language explanations.

## License

This project is open-source and available under the [MIT License](LICENSE).

---
*Note: 
- This system relies on the quality and coverage of the data loaded into Neo4j. If a reaction path exists in chemical literature but is not in the ORD subset loaded, the agent will not find it.
- If a target molecule cannot be synthesized using only the named starting materials, the system will not find a reaction path. Solvents, Catalysts and intermediate products do not have to be specified.*
