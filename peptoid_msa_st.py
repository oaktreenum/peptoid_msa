import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from io import BytesIO
import streamlit as st
import kaleido
import copy


# Allon Goldberg 
# Research Assistant, Flatiron Institute (NYC), Center for Computational Biology, Biomolecular Design Group
# 3/13/2025

# Description:
#     This script is a Streamlit-based web application designed for visualizing 
#     and analyzing multiple sequence alignments (MSA) of peptoids. The tool allows 
#     users to paste in sequences in FASTA format, apply custom color mappings to 
#     peptoid residues, and export visualizations in PNG or SVG formats.


def read_fasta(quasi_file):
    names = []
    seqs = []
    i_count = 0

    split_file = quasi_file.splitlines()
    for i in range(len(split_file)):
        if i == 0:
            if split_file[i].startswith(">"):  # Header line
                names.append(split_file[i][1:])
                i_count += 1
                seqs.append(split_file[i+1])
                i_count += 1
        if i >= i_count:
            if split_file[i].startswith(">"):  # Header line
                names.append(split_file[i][1:])
                i_count += 1
                seqs.append(split_file[i+1])
                i_count += 1
        else:
            continue
        
    return names, seqs

def plot_msa_plotly(names, seqs, color_map, property_map):
    seqs_array = np.array([[code for code in seq.split()] for seq in seqs])
    
    if seqs_array.ndim == 1:
        seqs_array = seqs_array.reshape(1, -1)

    len_seq, num_seq = seqs_array.shape
    ratio_len = 0.05
    ratio_num = 0.5
    w = num_seq * 43
    h = max(430,len_seq * 43)

    property_color_map = {}
    for idx, tuple in property_map.items():
        code = tuple[0]
        prop = tuple[1]
        if prop not in property_color_map:
            property_color_map[prop] = color_map.get(code, color_map['default'])

    # CREATE FIG
    fig = go.Figure()

    for i in range(seqs_array.shape[0]):
        for j in range(seqs_array.shape[1]):
            code = seqs_array[i, j] 
            color = color_map.get(code, color_map['default'])
            prop = property_map.get(code, 'default') 

            # Rectangles/boxes
            fig.add_trace(go.Scatter(
                x=[j, j + 1, j + 1, j],  
                y=[i, i, i + 1, i + 1],  
                fill='toself',            
                fillcolor=color,
                line=dict(color='rgba(0, 0, 0, 0)'),  
                mode='lines',             
                showlegend=False          
            ))

            # Text
            fig.add_trace(go.Scatter(
                x=[j + 0.5],             
                y=[i + 0.5],             
                text=[str(seqs_array[i, j])],  
                mode='text',                
                showlegend=False,
                textfont=dict(size=9, color='black')
            ))

    # Legend
    for prop, color in property_color_map.items():
        if prop != '—':
            fig.add_trace(go.Scatter(
                x=[None], y=[None],  
                mode='markers',
                marker=dict(size=10, color=color),
                name=prop 
            ))

    # Layout
    fig.update_layout(
        title='MSA Plot',
        margin=dict(l=134, r=88, t=100, b=88),
        xaxis=dict(
            tickmode='array',
            tickvals=np.arange(seqs_array.shape[1]) + 0.5,
            ticktext=[str(i + 1) for i in range(seqs_array.shape[1])],
            showgrid=False,            
            zeroline=False
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=np.arange(seqs_array.shape[0]) + 0.5,
            ticktext=names,
            showgrid=False,            
            zeroline=False,             
            autorange="reversed",
            tickfont=dict(size=13)      
        ),
        plot_bgcolor='white',         
        autosize=False,               
        width=max(300, ratio_num * num_seq * 888),  # Set plot width according to number of seqs
        height=min(0.3*max(300, ratio_num * num_seq * 888), max(300, ratio_len * len_seq * 818)),  # Set plot height according to length of seqs
    )

    return fig, w, h


# STREAMLIT SETUP
st.set_page_config(layout='wide', page_title = 'Peptoid MSA')
st.markdown('''
    <style>
        .block-container{
            width: 72%;  
        }
        .button-title{
            margin-bottom: 3px;
            font-weight: thin;
            color: white;
        }
        .center-button {
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100%;
        }
    </style>
''', unsafe_allow_html=True)

st.title('Peptoid MSA Visualization Tool')
st.markdown('<br>', unsafe_allow_html=True)
st.text('''This is a web-hosted python-based tool designed for visualizing multiple sequence alignments (MSAs) of peptoid sequences. It allows users to paste or upload peptoid sequences in FASTA format and generates a clear, interactive alignment plot. The tool supports customizable color schemes based on user-defined properties, helping to highlight key structural and chemical features.''')
st.markdown('''###### Click 'Show MSA' at bottom to view''')
st.markdown('<br>', unsafe_allow_html=True)


# TEXT INPUT FOR FASTA
example_fasta = '''> Sequence 1 \n208 PRO 208 632 314 208 632 332 624 632 631 601 HYP 634 HYP HYP 624 PRO 208 208 332 632 208 HYP HYP HYP 601 631 601 632 314 624 208 PRO 208 208 HYP 601 HYP 601 633 HYP 332 208 208 PRO 624 210 632
> Sequence 2 \n211 632 624 632 PRO 208 332 PRO 208 314 632 HYP HYP HYP 601 601 208 632 PRO 208 332 632 624 HYP HYP 633 HYP 601 601 208 332 624 208 PRO 208 202 HYP HYP HYP 632 601 601 PRO 208 208 332 624 211 632
> Sequence 3 \n208 632 624 PRO 332 208 632 PRO 208 332 HYP 601 HYP 631 HYP 601 208 332 632 208 632 PRO 624 HYP HYP 332 208 631 HYP 208 PRO 208 208 332 624 632 HYP HYP 631 HYP 208 601 PRO 624 208 332 203 632 632
> Sequence 4 \n211 332 624 632 314 208 PRO 632 624 PRO HYP 601 208 633 HYP HYP 624 632 632 208 332 PRO 210 HYP HYP HYP 601 631 601 632 332 208 211 PRO 624 211 HYP 631 HYP HYP 208 601 PRO 208 208 332 624 208 632'''
input_fasta = st.text_area('Paste FASTA sequence here', value=example_fasta, height=300)

st.markdown('<br>', unsafe_allow_html=True)

# CUSTOM CMAP MODULES
st.subheader('Define Custom Color Mapping')
# Define custom colormap
if 'custom_cmap' not in st.session_state:
    st.session_state.custom_cmap = {}
if 'custom_propmap' not in st.session_state:
    st.session_state.custom_propmap = {}
if 'standard_entries' not in st.session_state:
    st.session_state.standard_entries = {
        'chiral_hydrophobic': {'code': '601 602 621 622 623 624', 'prop': 'Chiral Hydrophobic', 'color': '#B95C00'},
        'hydrophobic': {'code': '001 003 005 007 020 101 103 127 130 202 203 208 210 211', 'prop': 'Hydrophobic', 'color': '#FFAF22'},
        'chiral_polar': {'code': '631 632 633 634', 'prop': 'Chiral Polar', 'color': '#0091B9'},
        'polar': {'code': '303 307 ', 'prop': 'Polar', 'color': '#88CFFF'},
        'polar_hydrophobic': {'code': '129', 'prop': 'Polar+Hydrophobic', 'color': '#C06EF7'},
        'negative': {'code': '314', 'prop': 'Negative', 'color': '#F95B5E'},
        'positive': {'code': '332 333', 'prop': 'Positive', 'color': '#5B75F9'},
        'pro': {'code': 'PRO', 'prop': 'Proline', 'color': '#B8B8B8'},
        'hyp': {'code': 'HYP', 'prop': 'Hydroxyproline', 'color': '#ABC5C5'},
        'sar': {'code': 'SAR', 'prop': 'Sarcosine', 'color': '#F9ECB3'},
        'default': {'code': 'default', 'prop': '', 'color': '#FFFFFF'}
    }
if 'skip_indices' not in st.session_state:
    st.session_state.skip_indices = []
if 'added_entries' not in st.session_state:
    st.session_state.added_entries = []
if 'rows' not in st.session_state:
    st.session_state.rows = {}
if 'rows_added' not in st.session_state:
    st.session_state.rows_added = {}

se_ids = ['chiral_hydrophobic', 'hydrophobic', 'chiral_polar', 'polar', 'polar_hydrophobic', 'negative', 'positive', 'pro', 'hyp', 'sar', 'default']
rows = st.session_state.rows
for i in range(len(se_ids)):
    col_obj = st.columns([2, 6, 2, 1])
    rows[i] = col_obj
st.session_state.rows = rows
se = st.session_state.standard_entries
for i, cols in rows.items():
    with cols[0]:
        se[se_ids[i]]['code'] = st.text_input('Enter 3-letter codes:', value=se[se_ids[i]]['code'], key=f'{se_ids[i]}_codes')
    with cols[1]:
        se[se_ids[i]]['prop'] = st.text_input('Describe properties for legend (ex: Chiral Hydrophobic):', value=se[se_ids[i]]['prop'], key=f'{se_ids[i]}_props')
    with cols[2]:
        se[se_ids[i]]['color'] = st.color_picker('Select color:', value=se[se_ids[i]]['color'], key=f'{se_ids[i]}_color')
st.session_state.standard_entries = se

for idx, entry in enumerate(st.session_state.standard_entries):
        code_vect = st.session_state.standard_entries[entry]['code'].split()
        for i in range(len(code_vect)):
            st.session_state.custom_cmap[code_vect[i]] = st.session_state.standard_entries[entry]['color']
            # st.session_state.custom_propmap[code_vect[i]] = st.session_state.standard_entries[entry]['prop']
            st.session_state.custom_propmap[idx] = [code_vect[i],st.session_state.standard_entries[entry]['prop']]   

rows_added = st.session_state.rows_added
for i in range(len(st.session_state.added_entries)-len(st.session_state.skip_indices)):
    col_obj = st.columns([2, 6, 2, 1])
    rows_added[i] = col_obj
st.session_state.rows_added = rows_added
i_row_added = 0
for idx, entry in enumerate(st.session_state.added_entries):
    if idx in st.session_state.skip_indices:
            continue
    with rows_added[i_row_added][0]:
        entry['code'] = st.text_input('Enter 3-letter codes:', key=f'code_input{idx}')
    with rows_added[i_row_added][1]:
        entry['prop'] = st.text_input('Describe properties for legend (ex: Chiral Hydrophobic):', key=f'prop_input{idx}')
    with rows_added[i_row_added][2]:
        entry['color'] = st.color_picker('Select color:', value=entry['color'], key=f'color_input{idx}') 
    with rows_added[i_row_added][3]:
        st.markdown("<p class='button-title'>Delete</p>", unsafe_allow_html=True)
        if st.button('x', key=f'delete_button{idx}'):
            st.session_state.skip_indices.append(idx)
            st.rerun()
    i_row_added += 1

if st.button('\+ Add Entry'):
    st.session_state.added_entries.append({
        'code': '',
        'prop': '',
        'color': '#ffffff'
    })
    st.rerun()

col1b, col2b, col3b, col4b = st.columns([2, 3, 1, 4])
with col1b:
    if st.button('Update Color Map'):
        for idx, entry in enumerate(st.session_state.standard_entries):
            code_vect = st.session_state.standard_entries[entry]['code'].split()
            for i in range(len(code_vect)):
                st.session_state.custom_cmap[code_vect[i]] = st.session_state.standard_entries[entry]['color']
                # st.session_state.custom_propmap[code_vect[i]] = st.session_state.standard_entries[entry]['prop']
                st.session_state.custom_propmap[idx] = [code_vect[i],st.session_state.standard_entries[entry]['prop']]
        for i in range(len(st.session_state.added_entries)):
            added_row = st.session_state.added_entries[i]
            st.session_state.custom_cmap[added_row['code']] = added_row['color']
            st.session_state.custom_propmap[i] = [added_row['code'], added_row['prop']]
with col4b:
    # Show the current color mappings
    if st.session_state.custom_cmap:
        st.subheader('Current Color Mapping:')
        with st.expander("Show Color Map", expanded=False):
            st.write(st.session_state.custom_cmap)

st.markdown('<br>', unsafe_allow_html=True)

def save_plot(fig, w, h, format="png"):
    buf = BytesIO()
    fig.write_image(buf, format=format, width=w, height=h, scale=3)
    return buf

st.markdown('# MSA Plot')

# MSA PLOT MODULES
if st.button('Show MSA'):
    if input_fasta[0] == '>':
        names, seqs = read_fasta(input_fasta)
        st.write('Analysing...')
        fig, w, h = plot_msa_plotly(names, seqs, st.session_state.custom_cmap, st.session_state.custom_propmap)
        print(w,h)
        st.plotly_chart(fig, use_container_width=False)
        col1c, col2c, col3c, col4c = st.columns(4)
        with col1c:
            st.write('Download buttons:')
        with col2c:
            st.download_button('PNG', save_plot(fig, w, h, 'png').getvalue(), 'msa.png', 'image/png', help='Click to download as .png')
        with col3c:
            st.download_button('SVG', save_plot(fig, w, h, 'svg').getvalue(), 'msa.svg', 'image/svg+xml', help='Click to download as .svg')
        with col4c:
            st.download_button('FASTA', data=input_fasta, file_name='sequences.fasta', mime='text/plain', help='Click to download sequences as fasta file')
    else:
        st.write('''Make sure you paste a FASTA sequence in the proper format, for example:\n
    > name1\n
    sequence1\n
    > name2\n
    seqence2\n
    > name3\n
    sequence3\n
    ...
    ''')


st.markdown('<br><br><br>', unsafe_allow_html=True)
    


