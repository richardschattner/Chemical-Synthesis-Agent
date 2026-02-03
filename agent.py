from langchain_core.runnables import RunnablePassthrough, RunnableLambda
from parser_agent import ChemicalParserAgent
from narrator_agent import ChemicalNarratorAgent
from retro import retrosynthesis_bfs, get_synthesis_subgraph, get_reaction_procedures, URI, USER, PASSWORD
from neo4j import GraphDatabase

# Initialize Neo4j Driver
driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))

def run_retrosynthesis(input_data):
    """
    Takes the output from ChemicalParserAgent and runs the BFS search.
    input_data: { "target": str, "starting_materials": List[str] }
    """
    if not input_data:
        return {"status": "error", "message": "Parsing failed, cannot run retrosynthesis."}

    target = input_data["target"]
    starting_materials = set(input_data["starting_materials"])
    
    print(f"Target: {target}")
    print(f"Building Blocks: \n{'\\n'.join(starting_materials)} provided")
    print(f"\n Starting Retrosynthesis Search...")
    # BFS
    path_rids = retrosynthesis_bfs(target, starting_materials, driver)

    if path_rids:
        # Get full details
        # Note: path_rids is [Last Step, ..., First Step] (Target -> Start)
        subgraph = get_synthesis_subgraph(path_rids, driver)
        return {
            "status": "found",
            "path_length": len(path_rids),
            "path_ids": path_rids, 
            "reactions": subgraph,
            "query": "Query logic encapsulated in get_synthesis_subgraph",
            "target": target,
            "starting_materials": list(starting_materials)
        }
    else:
        return {
            "status": "failed",
            "message": "No synthesis path found within max depth."
        }

def prepare_data_for_narrator(result):
    """
    If a path is found, retrieves the procedures in synthesized order
    and prepares the data for the narrator.
    """
    if result.get("status") != "found":
        return result

    # 1. Reverse path to get Synthesis Order (Start -> Target)
    # path_rids from BFS is [Target<-StepN, ..., Step1<-Start] (reverse synthesis)
    # We want [Step1, ..., StepN]
    synthesis_order_ids = result["path_ids"][::-1]
    
    # 2. Get procedures
    procedures = get_reaction_procedures(synthesis_order_ids, driver)
    
    # 3. Augment the procedures with structure info for the narrator
    # The 'reactions' list from get_synthesis_subgraph might be unordered or keyed differently.
    # We'll just pass the ordered procedures combined with the reaction details we already have.
    # Map reaction_id -> details
    reaction_map = {r['reaction']: r for r in result["reactions"]}
    
    annotated_steps = []
    for proc in procedures:
        rid = proc['reaction_id']
        details = reaction_map.get(rid, {})
        annotated_steps.append({
            "reaction_id": rid,
            "procedure_text": proc['procedure'],
            "reactants": details.get('reactants', []),
            "products": details.get('products', [])
        })
        
    return {**result, "annotated_steps": annotated_steps}

def format_structural_output(result):
    if result.get("status") == "error":
         return f"ERROR: {result.get('message')}"
         
    if result["status"] == "failed":
        return f"SEARCH FAILED: {result['message']}"

    output = [f"SUCCESS: Found path with {result['path_length']} reaction(s).\n"]
    
    output.append("--- Structural Data ---")
    for i, step in enumerate(result["reactions"], 1):
        output.append(f"Step {i}:")
        output.append(f"  Reaction ID: {step['reaction']}")
        output.append(f"  Reactants: {', '.join(step['reactants'])}")
        output.append(f"  Products: {', '.join(step['products'])}")
        output.append("")
    
    # Add Visualization Query
    if "path_ids" in result:
        output.append("--- Visualization Query (Neo4j) ---")
        output.append(f"MATCH (r:Reaction)-[rel]-(m:Molecule) WHERE r.reaction_id IN {str(result['path_ids'])} RETURN r, rel, m")
    
    return "\n".join(output)

# Setup the Chain
parser_agent = ChemicalParserAgent()
parser_runnable = parser_agent.get_runnable()

# Chain stops before narration/formatting string
retro_agent = (
    {"query": RunnablePassthrough()}
    | parser_runnable 
    | RunnableLambda(run_retrosynthesis)
    | RunnableLambda(prepare_data_for_narrator)
)

if __name__ == "__main__":
    print("Chemical Synthesis Agent")
    print("Enter your query:")
    user_query = input("Query: ").strip()

    if user_query:
        print(f"\nProcessing Query...")
        result = retro_agent.invoke(user_query)
        
        print("\n" + "="*30)
        # 1. Print Structural Info immediately
        print(format_structural_output(result))
        
        # 2. Stream Narrative if successful
        if result.get("status") == "found" and "annotated_steps" in result:
            print("\n=== AI GENERATED SYNTHESIS GUIDE ===")
            narrator = ChemicalNarratorAgent().get_runnable()
            # Stream the output
            for chunk in narrator.stream(result["annotated_steps"]):
                print(chunk, end="", flush=True)
            print("\n====================================\n")

    # Cleanup
    driver.close()
