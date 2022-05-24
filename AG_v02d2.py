# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 09:59:59 2021

@author: rmbranco
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
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
def gera_pop_aleat(i_min_val, i_max_val, i_size, p_size): 
    global individuos
    global fx
    individuos = []

    for i in range(p_size):
        crom_a = gera_cromossomo(i_min_val, i_max_val, i_size)
        individuos.append(crom_a)
        
    
def calcula_fx_pop(offset=0, faixa=5, decimais=2, n_bits=8):
    global l_fx
    global fx
    
    l_fx=[]
    fx=np.zeros(len(individuos))
    
    for i in range(len(individuos)):
        crom_dec = calc_decimal(individuos[i])
        y = calc_fx(crom_dec, offset, faixa, decimais, n_bits)
        l_fx.append(y)
    
    # converte fx para array
    fx = np.array(l_fx)

    
def calcula_fitness():
    global ftns
    global menor
    global maior
        
    maior = np.max(fx)    
    menor = np.min(fx)
    ftns = ((maior - fx) + menor) + 0.01
    
def calcula_prop():    
    global prop 
    prop = ftns / (np.sum(ftns))
      
def calcula_prop_ACC():    
    global prop_acc
    
    # criando array com proporcoes acumuladas
    p_size = len(individuos)
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
        
    novos_individuos = []
    
    for p in range(n_pares):
        xover(tmut, p)
   

def gera_novos_ind_aleat(i_min_val, i_max_val, i_size, num_indiv_aleat):
    for i in range(num_indiv_aleat):
        crom_a = gera_cromossomo(i_min_val, i_max_val, i_size)
        novos_individuos.append(crom_a)
        
def clonagem(num_indv_clones):
    global dict2
    
    dict1 = dict(zip(range(len(individuos)), fx))
    dict2 = sorted(dict1.items(), key=lambda x: x[1], reverse=False)
    for ind in range(num_indv_clones):
        pos_= int(dict2[ind][0])
        novos_individuos.append(individuos[pos_])

# funções gerais
def calcula_pesos(tamanho):
    global pesos
    pesos = np.zeros(tamanho, dtype=int)
    for p in range(tamanho):
        pesos[p]=np.power(2,p)

def mostra_pop():
    p_size = len(individuos)
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

def mostra_new_pop():
    # mostrando os dados produzidos 
    print("\n=======================================================================")
    print("= DADOS DA POPULAÇÃO GERADA === GEN =",num_gen, "==============================")
    print("=======================================================================")
    for j in range(len(individuos)):
        if j<2*num_pairs:
            print('Xover',j, individuos[j], '{0:.2f} {1:.3f} {2:.4f}'.format(fx[j], ftns[j], prop_acc[j]))
        if ((j>=2*num_pairs)and(j<(num_ind_aleat+2*num_pairs))):
            print('Aleat',j, individuos[j], '{0:.2f} {1:.3f} {2:.4f}'.format(fx[j], ftns[j], prop_acc[j]))
        if (j>=(num_ind_aleat+2*num_pairs)):
            print('Clone',j, individuos[j], '{0:.2f} {1:.3f} {2:.4f}'.format(fx[j], ftns[j], prop_acc[j]))
        
    print("=======================================================================\n")    

def copia_individuos():
    global individuos, novos_individuos
    
    individuos=[]
    individuos=novos_individuos
    novos_individuos=[]

# ----------------------------------------
# inicio do código
# ----------------------------------------

lote1={}
      
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
tam_pop = 30 # 30 indv no total
num_pairs = 9 # 18 individuos gerados por xover
taxa_mut = 1.0 # percentual = 1% no caso
num_clones = 2 # 20% do tam da população
num_ind_aleat = tam_pop - (2*num_pairs) - num_clones # restante é aleatorio
max_gen = 10 #num max de gerações
num_gen = 0 # num da geração corrente
evolut = []

# CHAMADA DAS FUNÇÕES PARA A PRIMEIRA POPULAÇÃO
calcula_pesos(tamanho)

# POP=00
gera_pop_aleat(l_min, l_max, tamanho, tam_pop)
calcula_fx_pop(offset, faixa, decimais, n_bits)

while True:
    
    calcula_fitness()
    calcula_prop()
    calcula_prop_ACC()
    best_ind = np.argmax(ftns)
    evolut.append(menor)
    print ("Gen:", num_gen, "/ Best indv: ", best_ind, "/ Best fx: ", fx[best_ind])    
    print ("Best solution: ", individuos[best_ind])
    h = np.histogram(fx, bins=10)
    print(h[0])
    
    if ((num_gen+1 == max_gen)or()):
        break

    sorteio_pais(tam_pop, num_pairs) 
    # mostra_selecionados(num_pairs)
    
    # xover c/ mutação
    gera_novos_ind_xover(num_pairs, 1.0)
    
    # cria novos ind aleatórios
    gera_novos_ind_aleat(l_min, l_max, tamanho, num_ind_aleat)
    
    # faz clonagem dos melhores
    clonagem(num_clones)
    
    copia_individuos()
    
    # calcula fx_individuos
    calcula_fx_pop(offset, faixa, decimais, n_bits)
      
    num_gen = num_gen + 1
        
fig, axs = plt.subplots(1, 2)
axs[0].plot(evolut)
x = [i for i in range(len(h[0]))]
axs[1].bar(x, h[0])
