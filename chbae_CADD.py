import streamlit as st
import pandas as pd
import numpy as np
import requests
from xml.etree import ElementTree as ET
from typing import List, Dict
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# PDB 그리기
from stmol import *
import py3Dmol
import time

#### 함수 단락
########################################
def get_reviewed_uniprot_ids_by_gene_symbols(gene_symbols):
    print("gene symbol을 통해 UniprotID를 반환합니다.")
    # """
    # gene symbols가 입력되면 
    # {"Gene symbol" : gene_symbol,
    #  "UniProt ID"  : uniprot_id,
    #  "UniProt sequence" : uniprot_sequence}
    # 를 원소로 가지는 list를 생성해서 df를 만들어 반환한다.
    # """    
    start_time = time.time()
    
    results = []
    
    for gene_symbol in gene_symbols:
        url = f"https://rest.uniprot.org/uniprotkb/search?query={gene_symbol}+AND+organism_id:9606+AND+reviewed:true&format=json"
        response = requests.get(url)
        
        if response.ok:
            data = response.json()
            found = False
            for protein_data in data.get('results', []):
                if any(gene.get('geneName', {}).get('value') == gene_symbol for gene in protein_data.get('genes', [])):
                    uniprot_id = protein_data['primaryAccession']
                    uniprot_seqeunce = protein_data['sequence']['value']
                    results.append({"Gene Symbol": gene_symbol, "UniProt ID": uniprot_id, "UniProt sequence": uniprot_seqeunce})
                    found = True
                    break  # 첫 번째 일치 항목 찾으면 종료
            if not found:
                # 일치하는 결과가 없을 때
                results.append({"Gene Symbol": gene_symbol, "UniProt ID": None, "UniProt sequence": uniprot_seqeunce})
        else:
            print(f"Failed to fetch data for {gene_symbol} from UniProt.")
            results.append({"Gene Symbol": gene_symbol, "UniProt ID": None, "UniProt sequence": uniprot_seqeunce})

        print(str(round(time.time() - start_time, 2)), "초가 걸렸습니다.")
    # DataFrame 생성
    df = pd.DataFrame(results)
    return df


def get_pdb_ids_from_uniprot_v2(uniprot_id):
    # print("UniprotID로 mutation 없고 wwPDB validation 결과를 반환합니다.")
    # """
    # uniprot_id를 RCSB PDB에 쿼리로 날려서 exact로 매치하면서 mutation이 없는 PDB ID리스트를 받아온다.
    # PDB에서 uniprot id로 수동 검색했을때 나오는 결과와 같다.
    # """
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": uniprot_id,
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": "UniProt",
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "equals",
                        "value": 0,
                        "attribute": "entity_poly.rcsb_mutation_count"
                    }
                }
            ]
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "polymer_entity"
    }

    response = requests.post(search_url, json=query)
    if response.status_code == 200:
        results = response.json()
        pdb_ids = [item["identifier"] for item in results["result_set"]]
        #print(len(pdb_ids), "개 검색되었습니다.")
        return pdb_ids if pdb_ids else None
    else:
        #print(uniprot_id)
        print("Failed to retrieve PDB IDs.")
        return None


def download_pdb(pdb_ids2):
    #print("PDB에서 구조 파일을 다운로드합니다.")
    # 검색 URL에 단순한 JSON 구조만 추가
    
    for id in pdb_ids2:
        print(id)
        pdb_url = f"https://files.rcsb.org/download/{id}.pdb"
        response = requests.get(pdb_url)
        if response.status_code == 200: 
            with open(f"./{id}.pdb", "wb") as file:
                file.write(response.content)
                print(f"Downloaded {id}.pdb")
        else:
            print(f"Failed to download {id}.pdb, error: {response.status_code}")


def extract_form_validation_xml(pdb_ids):
    print("wwPDB validation 파일을 확인합니다.")
    results = list()
    try:
        for index, id_ in enumerate(pdb_ids):
            # """
            # pdb_ids : 관심 Uniprot에 해당하는 pdb id가 담긴 list
            # validation report xml을 파싱해서 wwPDB 평가에 필요한 5개 항목을 반환함. 만약 값이 없으면 None으로 처리됨.            
            # """
            progress_bar.progress(int((index + 1) / len(pdb_ids) * 90), text='Progression rate')
            id_ = id_.lower()            
            id_split = id_.split("_")
            id_ = id_split[0]
            #print(id_, "처리중입니다.")
            validation_url = f"https://files.rcsb.org/pub/pdb/validation_reports/{id_[1:3]}/{id_}/{id_}_validation.xml.gz"
            tmp_req = requests.get(validation_url)
            #print(tmp_req.status_code)
            
            if hasattr(tmp_req, 'text'):
                xml_string = tmp_req.text
            else:
                xml_string = tmp_req
            
            root = ET.fromstring(xml_string)
            
            tag_name = 'Entry'
            attribute_names = ['absolute-percentile-DCC_Rfree',
                               'absolute-percentile-clashscore',
                               'absolute-percentile-percent-rama-outliers',
                               'absolute-percentile-percent-rota-outliers',
                               'absolute-percentile-percent-RSRZ-outliers',
                               'PDB-resolution',
                               'PDB-R',
                               'PDB-Rfree'                               
                              ]
            
            for element in root.findall(f'.//{tag_name}'):
                # 각 태그의 속성값들을 딕셔너리로 저장
                attributes = {}
                for attr_name in attribute_names:
                    # get() 메서드는 속성이 없는 경우 None을 반환
                    attributes[attr_name] = element.get(attr_name)
                results.append(attributes)
        return results
    except Exception as e:
        print(f'에러 발생: {str(e)}')


def extract_pdb_entity_seq(pdb_ids, uniprot_id):
    # """
    # pdb_ids: 6~7글자로 된 리스트. 함수 내에서 구분자로 쪼갠다음 url로 반환할 것임.
    # uniprot_id: 검색하려는 유전자의 대표 uniprot id
    # """
    try:
        tmp_list = list()
        pdb_ids_entity = [n.split("_") for n in pdb_ids]
        for entity, entity_id in pdb_ids_entity:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{entity}"
            req = requests.get(url)
            title_response = req.json()
            Entry_title = title_response['struct']['title']
            
            url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entity}/{entity_id}"
            #print(url)
            req = requests.get(url)
            entity_response  = req.json()
            #print(entity_response ['entity_poly']['rcsb_entity_polymer_type'])
            if entity_response['entity_poly']['rcsb_entity_polymer_type'] == "Protein": #DNA 나 다른 polymer일 수 있음
                if 'uniprot_ids' in entity_response['rcsb_polymer_entity_container_identifiers'].keys():
                    #print(entity_response['rcsb_polymer_entity_container_identifiers']['uniprot_ids'], uniprot_id)
                    if entity_response['rcsb_polymer_entity_container_identifiers']['uniprot_ids'][0] == uniprot_id:
                        #print("네", id_, entity_id, entity_response ['entity_poly']['pdbx_seq_one_letter_code_can'])
                        tmp_list.append(
                            {
                                'Entity' : entity,
                                'Entity_id' : entity_id,
                                'Entity_seq' : entity_response ['entity_poly']['pdbx_seq_one_letter_code_can'],
                                'Entry_title': Entry_title
                            }
                        )
                    else:
                        #print("아니오", entity, entity_id)
                        tmp_list.append(
                            {
                                'Entity' : entity,
                                'Entity_id' : entity_id,
                                'Entity_seq' : None,
                                'Entry_title': Entry_title
                            }
                        )

                #print("=======================")
            else:
                print(response.status_code)
                print("접속이 정상적이지 않습니다.")
        
        if tmp_list:
            df = pd.DataFrame(tmp_list)
            #print(df)
            return df
        else:
            print("없어요")
    except Exception as e:
        print(f'에러 발생: {str(e)}')

def where_is_sequence(full_seqeunce, want_to_know):
    # Perform local alignment
    alignments = pairwise2.align.localms(full_seqeunce, want_to_know, 2, -1, -0.5, -0.1)
    
    # Extract the best alignment
    best_alignment = alignments[0]
    
    # Print the alignment
    # print(format_alignment(*best_alignment))
    
    # Extract start and end indices for the alignment in the canonical sequence
    start_idx = best_alignment.start
    end_idx = best_alignment.end
    return start_idx, end_idx


def domain_extract_from_uniprot(uniprot_id, gene_symbol):
    # """
    # uniprot_id와 gene_symbol을 반환받아서 단백질의 도메인, 모티프 관련 정보를 따로 저장한다.
    # 'Uniprot ID': uniprot_id,
    # 'gene_symbol': gene_symbol,
    # 'type': type_,
    # 'start': location_start,
    # 'end': location_end,
    # 'description': description,
    # 'source': evidence.get('source'),
    # 'source_id': evidence.get('id')
    # """
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}'
    req = requests.get(url)
    response = req.json()

    result_list = list()
    #print(response)
    
    for item in response['features']:    
        if any(i in item['type'] for i in ['Region', 'Motif',
                                           'Site', 'DNA binding', 'Binding site']):
            type_ = item.get('type')
            location = item.get('location', {})
            location_start = location.get('start', {}).get('value')
            location_end = location.get('end', {}).get('value')
            description = item.get('description')
            ligand_name = item.get("ligand", {}).get("name")     
            if ligand_name:
                description = f"{description} (Ligand: {ligand_name})" if description else ligand_name
            evidences = item.get('evidences', [])
            for evidence in evidences:
                result_list.append(
                    {
                        'Uniprot ID': uniprot_id,
                        'gene_symbol': gene_symbol,
                        'type': type_,
                        'start': location_start,
                        'end': location_end,
                        'description': description,
                        'source': evidence.get('source'),
                        'source_id': evidence.get('id')
                    }
                )

    df = pd.DataFrame(result_list)
    return df


#### 사이드바에서 입력받아 데이터 처리하는 메인 함수
##############################
def main_function(gene_symbols):
    df = get_reviewed_uniprot_ids_by_gene_symbols(gene_symbols)
    uniprot_id = df['UniProt ID'].tolist()
    gene_name = df['Gene Symbol'].tolist()
    uniprot_seq = df['UniProt sequence'].tolist()

    summary = dict()
    summary2 = dict()
    for id_, name, canonical_sequence in zip(uniprot_id, gene_name, uniprot_seq):
        pdb_ids = get_pdb_ids_from_uniprot_v2(id_)
        if pdb_ids:
            pdb_ids_split = [n.split("_") for n in pdb_ids]
            pdb_ids2 = [n[0] for n in pdb_ids_split]
            print(pdb_ids2)
            #download_pdb(pdb_ids2) #다운 받는 라인
            results = extract_form_validation_xml(pdb_ids2) #Progress bar는 이 함수에서 작동함
            df2 = pd.DataFrame(results, index=pdb_ids)
            df2 = df2.apply(pd.to_numeric, errors='coerce')
            df2['average'] = df2[['absolute-percentile-DCC_Rfree',
                                  'absolute-percentile-clashscore',
                                  'absolute-percentile-percent-rama-outliers',
                                  'absolute-percentile-percent-rota-outliers',
                                  'absolute-percentile-percent-RSRZ-outliers']].mean(axis=1)
            df2['std'] = df2[['absolute-percentile-DCC_Rfree',
                              'absolute-percentile-clashscore',
                              'absolute-percentile-percent-rama-outliers',
                              'absolute-percentile-percent-rota-outliers',
                              'absolute-percentile-percent-RSRZ-outliers']].std(axis=1)
            df2 = df2[['absolute-percentile-DCC_Rfree',
                       'absolute-percentile-clashscore',
                       'absolute-percentile-percent-rama-outliers',
                       'absolute-percentile-percent-rota-outliers',
                       'absolute-percentile-percent-RSRZ-outliers',
                       'average',
                       'std',
                       'PDB-resolution',
                       'PDB-R',
                       'PDB-Rfree']]
                    
            # """
            # 여기에 PDB ID 별 시퀀스 비교해서 어떤 도메인인지 확인하는 기능이 필요하네. 그래서 df로 취합해야 엑셀로 정리 가능.
            # """
            df3 = extract_pdb_entity_seq(pdb_ids, id_)
            df2['Entity'] = df3['Entity'].values #인덱스가 안맞아서 강제로 할당함
            df2['Entity_id'] = df3['Entity_id'].values
            df2['Entity_seq'] = df3['Entity_seq'].values
            df2['Title'] = df3['Entry_title'].values
    
            #아래부터는 uniprot_seq와 각 entity의 seq 위치를 비교해서 position을 저장하는 단락. 함수로 안빼고 그냥 바로 처리함.
            position = list()
            for seq in df2['Entity_seq']:
                if pd.isna(seq):  # 빈 값 처리
                    position.append(np.nan)
                else:
                    start_index, end_index = where_is_sequence(canonical_sequence, seq)
                    position.append(f"{start_index}-{end_index}")
            df2['position'] = position
            print(df2['position'], position)
    
            df2 = df2[['absolute-percentile-DCC_Rfree',
                       'absolute-percentile-clashscore',
                       'absolute-percentile-percent-rama-outliers',
                       'absolute-percentile-percent-rota-outliers',
                       'absolute-percentile-percent-RSRZ-outliers',
                       'average',
                       'std',
                       'PDB-resolution',
                       'PDB-R',
                       'PDB-Rfree',
                       'Entity',
                       'Entity_id',
                       'Entity_seq',
                       'position',
                       'Title']]
    
            #summary[f"{id_}_{name}_{str(len(canonical_sequence))}"] = df2
            summary[name] = [id_, canonical_sequence, df2]
        else:
            print("없어요")
    
       
        #v1.4부터 추가한 단백질별 모티프, 도메인 정보 따로 저장하기. 위 데이터와 병합하기에는 구조가 너무 다르다.
        df4 = domain_extract_from_uniprot(id_, name)
        #summary2[f"{id_}_{name}_{str(len(canonical_sequence))}"] = df4
        summary2[name] = [id_, canonical_sequence, df4]
        
    return summary, summary2
############################################################
############################################################

##### Page layout setting
##############################
st.set_page_config(
    page_title="Protein Structure informatics",
    layout="wide"
)   

if "initialized" not in st.session_state:
    st.session_state["initialized"] = False
    st.session_state["form_data"] = None

##### 사이드바
##############################
with st.sidebar:
    st.write('ch Bae')
    st.write("""
    탐색할 유전자심볼을 입력해주세요. 여러개를 입력할땐 줄 나눔을 해주세요.
    """)

    with st.form('form1'):
        txt = st.text_area(
            'Input gene symbol list',
            value = '',
            height=500,
        )
        gene_symbols = txt.split('\n')
        submitted = st.form_submit_button('submit')
        if submitted:
            st.session_state["initialized"] = True
            st.session_state['gene_symbols'] = gene_symbols
            st.write('Submitted.')


##### 메인 단락
########################################
st.title('Protein Structure Informatics')
st.subheader('This application provides a feature that quickly summarizes protein structures using Gene symbols, UniProt, and PDB.')

if (st.session_state["initialized"] == True) and ('gene_symbols' in st.session_state):
    progress_bar = st.progress(0, text='Progression rate')
    summary, summary2 = main_function(gene_symbols)
    gene_symbol = st.selectbox(
        "원하는 gene symbol에 대한 결과를 보입니다",
        (gene_symbols),
    )
    st.write(f'Canonical protein sequence of {gene_symbol} (UniProt ID {summary[gene_symbol][0]}) is below ({len(summary[gene_symbol][0])} amino acids):')
    canonical_sequence = summary[gene_symbol][1]
    canonical_sequence_display = ' '.join([canonical_sequence[i:i+5] for i in range(0, len(canonical_sequence), 5)])
    wrapped = '\n'.join([canonical_sequence_display[i:i+60] for i in range(0, len(canonical_sequence_display), 60)])
    st.code(wrapped, line_numbers=True)
   
    st.write(f'{gene_symbol}의 PDB accession table')
    st.write('Index is represented as "PDB ID_Entity ID"')
    st.dataframe(summary[gene_symbol][2], use_container_width=True)
    st.write(f'{gene_symbol}의 Structure info table')
    st.dataframe(summary2[gene_symbol][2], use_container_width=True)

    progress_bar.progress(100, 'Done')#완료가 됐다는 느낌을 주기 위해

    ##### py3Dmol단락
    pdb_id_tmp = list(summary[gene_symbol][2].index)
    pdb_to_render = [item.split('_') for item in pdb_id_tmp]
    pdb_to_render = [n[0] for n in pdb_to_render]
    pdb_display = st.selectbox(
        "원하는 gene symbol에 대한 결과를 보입니다",
        (pdb_to_render),
        index=None,
        placeholder="Please wait until PDB is loading."
    )
    if pdb_display != None:
        st.write('* 선택한 PDB ID :', pdb_display, '마우스 휠로 확대/축소 가능하고 드래그로 방향을 바꿀 수 있습니다.')
        xyzview = py3Dmol.view(query='pdb:'+pdb_display) 
        xyzview.setStyle({'cartoon':{'color':'spectrum'}})
        showmol(xyzview, height =500, width=1000)
