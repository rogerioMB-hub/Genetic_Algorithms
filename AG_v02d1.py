# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 09:59:59 2021

@author: rmbranco
"""

import numpy as np

# funções do individuo      
def gera_cromossomo(min_, max_, size_):
    return np.random.randint(min_, max_+1, size_)    

def calc_decimal(crom_):
    return np.dot(crom_, pesos)        

def calc_fx(valor_d, offset_, faixa, n_casas, n_bits):
    aux1 = np.round(faixa*(10**n_casas))
    resol = (2**n_bits)-1           
    x = offset_ + ((aux1 * valor_d / resol)/(10.0**n_casas))
    return x**2

# funções da população
def gera_pop_aleat(i_min_val, i_max_val, i_size, p_size, offset=0, faixa=5, decimais=2, n_bits=8):
    global individuos
    global fx
    individuos = []

    l_fx=[]
    for i in range(p_size):
        crom_a = gera_cromossomo(i_min_val, i_max_val, i_size)
        individuos.append(crom_a)
        crom_dec = calc_decimal(crom_a)
        y = calc_fx(crom_dec, offset, faixa, decimais, n_bits)
        l_fx.append(y)
    # converte fx para array
    fx = np.array(l_fx)
          
def calcula_fitness():
    global ftns
        
    maior = np.max(fx)    
    menor = np.min(fx)
    ftns = ((maior - fx) + menor) + 0.01
    
def calcula_prop():    
    global prop 
    prop = ftns / (np.sum(ftns))
      
def calcula_prop_ACC(p_size):    
    global prop_acc
    
    # criando array com proporcoes acumuladas
    prop_acc = np.zeros(p_size)
      
    prop_acc[0] = prop[0]
    for i in range(1,p_size):
        prop_acc[i] = prop_acc[i-1]+prop[i]
    prop_acc = prop_acc * 360.0
    
def sorteio_pais(p_size, num_pairs):
    global pais
    pais = np.zeros((num_pairs,2), dtype=int)
  
    for s in range(num_pairs):
        while pais[s,0]==pais[s,1]:
            lances = np.random.random_sample((2,)) * 360.0
            for p in range(2):
               for k in range(p_size):
                  if prop_acc[k]>lances[p]:
                     pais[s,p]=k
                     break

def calcula_pto_corte(size_crom):
    return np.random.randint(1, size_crom)

def flip(valor):
    if valor==0:
        valor=1
    else:   
        valor=0
    
    return valor

def mask_mutation(size_crom):
    return np.random.random_sample(size_crom)*100

def xover(tmut, pair):
    p1 = individuos[pais[pair,0]]
    p2 = individuos[pais[pair,1]]   
    
    size_crom = np.shape(p1)[0]
    corte = calcula_pto_corte(size_crom)
    
    f1 = np.zeros(size_crom, dtype=int)
    f2 = np.zeros(size_crom, dtype=int)
    
    # f1=np.concatenate([p1[:corte], p2[corte:]])
    # f2=np.concatenate([p2[:corte], p1[corte:]])

    mascara=mask_mutation(size_crom)
    for i in range(size_crom):
        if (i<corte):
            f1[i]=p1[i]
            f2[i]=p2[i]
        else:
            f1[i]=p2[i]
            f2[i]=p1[i]
        # etapa de mutação 
        if (mascara[i]<tmut):
            f1[i]=flip(f1[i])
            f2[i]=flip(f2[i])
    
    novos_individuos.append(f1)
    novos_individuos.append(f2)

def gera_novos_ind_xover(n_pares, tmut):
    global novos_individuos
    global nfx
    global nftns
    
    novos_individuos = []
    
    for p in range(n_pares):
        xover(tmut, p)
    
    nl_fx=[]
    for i in range(len(novos_individuos)):
        crom_dec = calc_decimal(novos_individuos[i])
        y = calc_fx(crom_dec, offset, faixa, decimais, n_bits)
        nl_fx.append(y)
    # converte fx para array
    nfx = np.array(nl_fx)    
    maior = np.max(nfx)    
    menor = np.min(nfx)
    nftns = ((maior - nfx) + menor) + 0.01
    
             
# funções gerais
def calcula_pesos(tamanho):
    global pesos
    pesos = np.zeros(tamanho, dtype=int)
    for p in range(tamanho):
        pesos[p]=np.power(2,p)

def mostra_pop(p_size):
    # mostrando os dados produzidos 
    print("\n=======================================================================")
    print("= DADOS DA POPULAÇÃO GERADA ===========================================")
    print("=======================================================================")
    for j in range(p_size):
        print(j, individuos[j], '{0:.2f} {1:.3f} {2:.4f}'.format(fx[j], ftns[j], prop_acc[j]))
        #print(j, individuos[j], fx[j], ftns[j], prop_acc[j])
    print("=======================================================================\n")
    
def mostra_selecionados(n_pares):   
    print("\n=======================================================================")
    print("= INDIVIDUOS SELECIONADOS PARA CRUZAMENTO =============================")
    print("=======================================================================")
    
    for k in range(n_pares):
        print("Pais ", k, pais[k])
    print("=======================================================================")    

def mostra_filhos():
    # mostrando os dados produzidos 
    print("\n=======================================================================")
    print("= DADOS DA POPULAÇÃO GERADA POR xover ===================================")
    print("=======================================================================")
    for j in range(len(novos_individuos)):
        print(j, novos_individuos[j], '{0:.2f} {1:.3f}'.format(nfx[j], nftns[j]))
        
    print("=======================================================================\n")    
# ----------------------------------------
# inicio do código
# ----------------------------------------
      
# VARIÁVEIS 
#    do individuo
l_min = 0
l_max = 1
tamanho = 8

#    da função
offset=0
faixa=5.0
decimais=2
n_bits=8

#    da população
tam_pop = 10
num_pairs = 3
taxa_mut = 1.0 # percentual = 1% no caso

# CHAMADA DAS FUNÇÕES PARA A PRIMEIRA POPULAÇÃO
calcula_pesos(tamanho)

# POP=00
gera_pop_aleat(l_min, l_max, tamanho, tam_pop, offset, faixa, decimais, n_bits)
calcula_fitness()
calcula_prop()
calcula_prop_ACC(tam_pop)

mostra_pop(tam_pop)

sorteio_pais(tam_pop, num_pairs) 
mostra_selecionados(num_pairs)

# xover c/ mutação
gera_novos_ind_xover(num_pairs, 1.0)
mostra_filhos()
# composiçào da nova população

# se obj n atendidos, retoma calculo de fitness e demais passos
