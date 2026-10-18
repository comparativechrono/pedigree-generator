import streamlit as st
from PedGen import pedigree_generator, construct_pedigree_graph, plot_pedigree_tree

st.title('Pedigree Generator')


#sliders, buttons, and select boxes
st.text_input("Family ID", key='FamilyID')
max_children = st.slider(label= 'Max Children',
                         min_value= 1,
                         max_value= 10,
                         value= 3,
                         help= 'The maximum number of children each spousal pair can have')

mode = st.selectbox(label= 'Mode of Inheritance',
                    options= ['AD', 'AR'],
                    index= 0,
                    format_func= lambda x: 'Autosomal Dominant' if x == 'AD' else 'Autosomal Recessive',
                    help= 'The mode of inheritence dictating the phenotype inheritance pattern')

generation_count = st.slider(label= 'Number of Generations',
                             min_value= 1,
                             max_value= 10,
                             value= 3)
AffectedSpouse = st.checkbox(label = 'Affected Secondary Founders',
                             help= 'Whether the partners who marry into the family can carry the alternate allele')
if AffectedSpouse:
    alt_freq = st.slider(label= 'Alternate Allele Frequency',
                         min_value=0.01,
                         max_value=1.00,
                         step=0.01,
                         value= 0.10,
                         help='How common the alternate allele is in the population that partners come from')
else:
    alt_freq = 0

AdvancedSettings = st.toggle(label = 'Advanced Settings',
                             help= 'Set how often partners have parents in the pedigree, and how often children have partners')
if AdvancedSettings:
    BackpropLikelihood = st.slider(label= 'Ancestry Backprop Likelihood',
                                   min_value= 0.0,
                                   max_value=1.0,
                                   step=0.1,
                                   value=0.5)
    SpouseLikelihood = st.slider(label= 'Reproductive Likelihood',
                                 min_value= 0.0,
                                 max_value=1.0,
                                 step=0.1,
                                 value=0.5)
else:
    BackpropLikelihood = 0.3
    SpouseLikelihood = 0.6

FamilyID = st.session_state.FamilyID

generated_df = pedigree_generator(max_children= max_children,
                             FamilyID= FamilyID,
                             mode= mode,
                             generation_count= generation_count,
                             AffectedSpouse= AffectedSpouse,
                             alt_freq= alt_freq,
                             BackpropLikelihood= BackpropLikelihood,
                             SpouseLikelihood= SpouseLikelihood)

st.write(plot_pedigree_tree(generated_df, f'{FamilyID} Pedigree' if FamilyID else 'Randomly Generated Pedigree'))
st.dataframe(
    generated_df,
    hide_index= True,
    column_order= ['FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'],
    column_config= {
        'FamilyID': 'Family ID', 
        'IndividualID': 'ID', 
        'PaternalID': 'Paternal ID', 
        'MaternalID': 'Maternal ID',
    }
)

pedfile = generated_df.to_csv(columns= ['FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'],
                              header= False,
                              index= False,
                              sep= ' ')
st.download_button(
    label= "Download PedFile",
    data= pedfile,
    file_name=f'{FamilyID}.ped' if FamilyID else 'GeneratedFamily.ped',
    mime= "text/plain",
    icon=':material/download:')