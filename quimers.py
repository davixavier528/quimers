#!/usr/bin/env python3

"""
utilização:
python3 quimers.py adjuvants_list.fasta mhci_epitopes_list linker_mhci mhcii_epitopes_list linker_mhcii his_tag

o arquivo adjuvants_list.fasta precisa estar no formato fasta linear
o arquivo mhci_epitopes_list precisa ser unicamente composto pelo epítopos MHCI
o arquivo mhcii_epitopes_list precisa ser unicamente composto pelo epítopos MHCII
o argumento his_tag precisa ser 'yes' ou 'no'
"""

# importar bibliotecas necessárias
import sys
import random
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# criar variáveis de acordo com o input do usuário
adjuvant       = sys.argv[1]
mhci_epitopes  = sys.argv[2]
mhci_linker    = sys.argv[3]
mhcii_epitopes = sys.argv[4]
mhcii_linker   = sys.argv[5]
his_tag        = sys.argv[6]

# linker do adjuvante
adj_linker = 'EAAAK'

# definir localização dos epítopos mais próximos da região N-terminal
# prior = 'mhci', os epítopos MHCI são os primeiros
# prior = 'mhcii', os epítopos MHCII são os primeiros
prior = 'mhci'

# quantidade de iterações do Monte Carlo
max_iter = 10000

# definir função que calcula índice alifático
# baseado em Ikai (1980)
def aliphatic_index(protein):
    alip_indx = "%0.2f" % ((protein.get_amino_acids_percent()['A']*100)+(2.9*(protein.get_amino_acids_percent()['V']*100))+3.9*((protein.get_amino_acids_percent()['I'])*100+(protein.get_amino_acids_percent()['L'])*100))
    return alip_indx

# ler arquivos de texto
with open(adjuvant,'r') as a, open(mhci_epitopes, 'r') as b, open(mhcii_epitopes, 'r') as c:

    adj_list       = a.read().splitlines()
    mhci_epi_list  = b.read().splitlines()
    mhcii_epi_list = c.read().splitlines()

    # processar adjuvantes
    for line in adj_list:
        if line.startswith('>'):
            adj_name = line.replace('>','')
        else:
            # criar lista vazia que armazena os menores valores de ii
            ii_log = [ ]

            # criar lista vazia que armazena as sequências das proteínas quimeras que apresentaram os menores valores de ii
            quim_prot_log = [ ]

            # iniciar Monte Carlo
            for i in range(max_iter):

                # gerar compilação aleatória dos epítopos e proteína quimérica
                if mhci_linker == mhcii_linker:
                    epitopes = mhci_epi_list + mhcii_epi_list
                    random.shuffle(epitopes)
                    quim_epitopes = mhci_linker.join(epitopes)
                    quim_prot = line + adj_linker + quim_epitopes
                elif mhci_linker != mhcii_linker:
                    random.shuffle(mhci_epi_list)
                    quim_mhci  = mhci_linker.join(mhci_epi_list)
                    random.shuffle(mhcii_epi_list)
                    quim_mhcii = mhcii_linker.join(mhcii_epi_list)
                    if prior == 'mhci':
                        quim_prot = line + adj_linker + quim_mhci + mhcii_linker + quim_mhcii
                    elif prior == 'mhcii':
                        quim_prot = line + adj_linker + quim_mhcii + mhci_linker + quim_mhci
                
                # adicionar his-tag
                if his_tag == 'yes':
                    quim_prot += 'HHHHHH'

                # calcular o índice de instabilidade
                quim_prot_analysis = ProteinAnalysis(quim_prot)
                ii = "%0.4f" % quim_prot_analysis.instability_index()

                # armazenar o primeio resultado do for loop
                if len(ii_log) == 0 and len(quim_prot_log) == 0:
                    ii_log.append(ii)
                    quim_prot_log.append(quim_prot)

                # armazenar somente os menores índices de instabilidade
                elif ii < ii_log[-1]:
                    del ii_log[-1]
                    del quim_prot_log[-1]
                    ii_log.append(ii)
                    quim_prot_log.append(quim_prot)
            
            # calcular propriedades da melhor conformação para o adjuvante
            final_quim_analysis = ProteinAnalysis(quim_prot_log[-1])
            qtd_aa = len(quim_prot_log[-1])
            mw = "%0.2f" % final_quim_analysis.molecular_weight()
            ip = "%0.2f" % final_quim_analysis.isoelectric_point()
            final_ii = ii_log[-1]
            ai = aliphatic_index(final_quim_analysis)
            gravy = "%0.4f" % final_quim_analysis.gravy()

            # exibir resultados finais
            print ( )
            print ( '# A conformação com o melhor índice de instabilidade encontrada para o adjuvante "' + adj_name + '" possui as seguintes características:' )
            print ( '# Quantidade de aminoácidos:' , qtd_aa )
            print ( '# Peso molecular:' , mw )
            print ( '# Ponto isoelétrico:' , ip )
            print ( '# Índice de instabilidade:' , final_ii )
            print ( '# Índice alifático:' , ai )
            print ( '# GRAVY:' , gravy )
            print ( '# Sequência:' )
            print ( '>quimeric_protein_' + adj_name + '_adjuvant' )
            print ( quim_prot_log[-1] )