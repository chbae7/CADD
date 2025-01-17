import streamlit as st
import pandas as pd
import numpy as np


##### 사이드바
##############################
with st.sidebar:
    st.write("""
    탐색할 유전자심볼을 입력해주세요. 여러개를 입력할땐 줄 나눔을 해주세요.
    """)

    with st.form('form1'):
        txt = st.text_area(
            'Input gene symbol list',
            value = '',
            height=500,
        )
        txt_to_list = txt.split('\n')
        st.write(txt_to_list)

        submitted = st.form_submit_button('submit')
        if submitted:
            st.write('섭밋됐습니다.')

##### 메인 단락
########################################
st.write('HELLO WORLD')
st.write('RESULT')
option = st.selectbox(
    "원하는 gene symbol에 대한 결과를 보입니다",
    (txt_to_list),
)

st.write('CASP_A00000_234 과 같은 식으로 출력해야함')

canonical_sequence = 'MQMSPALTCLVLGLALVFGEGSAVHHPPSYVAHLASDFGVRVFQQVAQASKDRNVVFSPYGVASVLAMLQLTTGGETQQQIQAAMGFKIDDKGMAPALRHLYKELMGPWNKDEISTTDAIFVQRDLKLVQGFMPHFFRLFRSTVKQVDFSEVERARFIINDWVKTHTKGMISNLLGKGAVDQLTRLVLVNALYFNGQWKTPFPDSSTHRRLFHKSDGSTVSVPMMAQTNKFNYTEFTTPDGHYYDILELPYHGDTLSMFIAAPYEKEVPLSALTNILSAQLISHWKGNMTRLPRLLVLPKFSLETEVDLRKPLENLGMTDMFRQFQADFTSLSDQEPLHVAQALQKVKIEVNESGTVASSSTAVIVSARMAPEEIIMDRPFLFVVRHNPTGTVLFMGQVMEP'
canonical_sequence_display = ' '.join([canonical_sequence[i:i+5] for i in range(0, len(canonical_sequence), 5)])
wrapped = '\n'.join([canonical_sequence_display[i:i+60] for i in range(0, len(canonical_sequence_display), 60)])
st.code(wrapped, line_numbers=True)

## 데이터프레임 출력 단락
st.write('CASP3의 PDB accession table1')
df1 = pd.read_excel('wwPDBvalidation_테스트.xlsx')
st.dataframe(df1)

st.write('CASP3의 PDB accession table2  이건 도메인 같은 정보임')
df2 = pd.read_excel('wwPDBvalidation_테스트_structure_info.xlsx')
st.dataframe(df2)
