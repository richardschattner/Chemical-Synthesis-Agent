import streamlit as st
import networkx as nx
from pyvis.network import Network
import tempfile
import os
import json

# Import the backend logic
from agent import retro_agent
from narrator_agent import ChemicalNarratorAgent

st.set_page_config(page_title="Chemical Synthesis Agent", layout="wide")

st.title("🧪 Chemical Synthesis AI")

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Initialize last graph
if "graph_html" not in st.session_state:
    st.session_state.graph_html = None

def generate_graph(reactions, target=None, starting_materials=None):
    """
    Generates an interactive PyVis graph from the reaction list.
    """
    if starting_materials is None: starting_materials = []
    
    # Transparent background (using #000000 essentially or rgba(0,0,0,0))
    # We use a black background to look good in dark mode if transparency fails in iframe behavior
    net = Network(height="600px", width="100%", bgcolor="#222222", font_color="white", directed=True)
    
    # Global Options for styling
    options = {
        "nodes": {
            "font": {
                "size": 16, # Larger font -> Larger nodes
                "color": "white",
                "strokeWidth": 0, # No black border on text
                "face": "arial"
            },
            "borderWidth": 2
        },
        "edges": {
            "color": {"inherit": True},
            "smooth": False,
            "font": {
                "size": 12,
                "color": "white",
                "strokeWidth": 0,
                "align": "horizontal"
            }
        },
        "physics": {
            "forceAtlas2Based": {
                "gravitationalConstant": -100,
                "centralGravity": 0.01,
                "springLength": 200,
                "springConstant": 0.08
            },
            "minVelocity": 0.75,
            "solver": "forceAtlas2Based"
        }
    }
    net.set_options(json.dumps(options))
    
    # Track added nodes to avoid duplicates
    added_nodes = set()
    
    for step in reactions:
        rid = step['reaction']
        
        # Add Reaction Node (Box puts label inside)
        if rid not in added_nodes:
            # Shorten Reaction ID for display
            display_rid = rid[:10] + "..." if len(rid) > 10 else rid
            # Use 'box' to enforce label inside
            net.add_node(rid, label=display_rid, color="#ff7f50", shape="box", title=f"Reaction: {rid}")
            added_nodes.add(rid)
            
        # Add Reactants -> Reaction
        for r in step['reactants']:
            if r not in added_nodes:
                # Color Logic
                color = "#4ECDC4" # Default Cyan
                if r in starting_materials:
                    color = "#1E90FF" # Dodger Blue
                elif r == target:
                    color = "#32CD32" # Lime Green
                
                # Truncate label to keep circle distinct and not too wide
                label = r[:10] + "..." if len(r) > 10 else r
                # Use 'circle' to enforce label inside
                net.add_node(r, label=label, title=r, color=color, shape="circle")
                added_nodes.add(r)
            
            # Edge: Molecule -> Reaction
            net.add_edge(r, rid, label="Reactant", color="#aaaaaa")
            
        # Reaction -> Products
        for p in step['products']:
            if p not in added_nodes:
                # Color Logic
                color = "#4ECDC4" 
                if p == target:
                    color = "#32CD32" 
                elif p in starting_materials:
                    color = "#1E90FF" 
                
                label = p[:10] + "..." if len(p) > 10 else p
                # Use 'circle'
                net.add_node(p, label=label, title=p, color=color, shape="circle")
                added_nodes.add(p)
                
            # Edge: Reaction -> Molecule
            net.add_edge(rid, p, label="Product", color="#aaaaaa")
            
    # Save to temp file and read back
    try:
        fd, path = tempfile.mkstemp(suffix=".html")
        os.close(fd)
        net.save_graph(path)
        with open(path, 'r', encoding='utf-8') as f:
            html = f.read()
            
        # INJECT CSS TO FORCE BACKGROUND COLOR #222222
        # PyVis sometimes only sets the network div background, leaving body white.
        # We also remove the border
        html = html.replace('</head>', '<style>body { background-color: #222222 !important; margin: 0; padding: 0; } #mynetwork { border: none !important; }</style></head>')
            
        os.remove(path)
        return html
    except Exception as e:
        return f"Error generating graph: {e}"

# Layout: Split screen
col1, col2 = st.columns(2)

with col1:
    st.subheader("Assistant Chat")
    
    # render history
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    # Input
    if prompt := st.chat_input("Enter synthesis target (e.g. 'Make aspirin from ...')"):
        
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)

        with st.chat_message("assistant"):
            message_placeholder = st.empty()
            full_response = ""
            
            with st.spinner("Analyzing Chemical Space..."):
                try:
                    # Invoke Backend
                    result = retro_agent.invoke(prompt)
                    
                    if result.get("status") == "found":
                        # Prepare Narrative
                        narrator = ChemicalNarratorAgent().get_runnable()
                        
                        # Generate Graph
                        graph_html = generate_graph(
                            result["reactions"], 
                            target=result.get("target"), 
                            starting_materials=result.get("starting_materials")
                        )
                        st.session_state.graph_html = graph_html
                        
                        # Stream Narrative
                        if "annotated_steps" in result:
                            stream = narrator.stream(result["annotated_steps"])
                            for chunk in stream:
                                full_response += chunk
                                message_placeholder.markdown(full_response + "▌")
                            message_placeholder.markdown(full_response)
                        else:
                            full_response = "Found path, but unable to generate narrative."
                            message_placeholder.markdown(full_response)
                            
                    else:
                        full_response = f"❌ Search Failed: {result.get('message', 'No path found.')}"
                        message_placeholder.markdown(full_response)
                        # Don't clear graph on failure, maybe user wants to see previous? 
                        # Or maybe clear it? Let's keep it for context unless successful new one replaces it.
                
                except Exception as e:
                    full_response = f"⚠️ An error occurred: {str(e)}"
                    message_placeholder.markdown(full_response)
                    
            st.session_state.messages.append({"role": "assistant", "content": full_response})
            # Force rerun to update the graph on the right immediately (sometimes needed)
            st.rerun()

with col2:
    st.subheader("Synthesis Graph")
    
    # Legend
    st.markdown("""
    <style>
        .legend-box { display: inline-block; width: 15px; height: 15px; margin-right: 5px; vertical-align: middle; }
    </style>
    <div style="margin-bottom: 10px; background-color: #222; padding: 10px; border-radius: 5px; color: #ddd;">
        <span class="legend-box" style="background-color: #1E90FF;"></span>Starting Molecule &nbsp;&nbsp;
        <span class="legend-box" style="background-color: #32CD32;"></span>Target Product &nbsp;&nbsp;
        <span class="legend-box" style="background-color: #4ECDC4;"></span>Intermediate/Other &nbsp;&nbsp;
        <span class="legend-box" style="background-color: #ff7f50;"></span>Reaction
    </div>
    """, unsafe_allow_html=True)

    if st.session_state.graph_html:
        st.components.v1.html(st.session_state.graph_html, height=750, scrolling=True)
    else:
        st.info("Enter a query to visualize the reaction pathway.")
        st.warning("""⚠️ **Note:** Synthesis paths of the target molecule are found based on the provided starting components.
            If the starting components are not sufficient to synthesize the target, the search will fail. 
            Solvents and catalysts are inferred from the database and do not need to be mentioned explicitly in your query.""")
        st.markdown("""**Try queries like:** \n
        - *How can I synthesize aspirin from salicylic acid and acetic anhydride?*\n
        - *Synthesize paracetamol from 4-aminophenol and acetic anhydride.*\n
        - *Find a path to make Benzocaine using p-aminobenzoic acid and ethanol.*""")