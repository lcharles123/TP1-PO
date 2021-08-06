import sys
from fractions import Fraction
DEBUG1=False # tabela da PL aux pronta pro simplex 
DEBUG2=False # tabela da PL aux apos o simplex  
DEBUG3=False # tabela da PL principal antes do simplex
DEBUG4=True # tabela da PL principal apos do simplex
DEBUG5=False # tabela da PL principal ilimitada

'''Transpoe para uma matriz representada por lista de listas'''
def transpor(M): 
    return list(map(list, zip(*M)))

'''Recebe uma coluna e retorna o indice da base, ou -1 caso nao seja'''
def base(v): 
    indice = -1 
    for i in range(len(v)):
        if v[i] == 0:
            continue
        elif v[i] == 1 and indice == -1:
            indice = i 
        else:
            return -1
    return indice

'''Monta a PL auxiliar e executa o simplex, retorna o tableau resultante com uma base viavel, caso otima ou ilimitada'''
def fase2(n, m, d, cAUX, M, A, b):
    fobj = 0
    A = transpor(A)             # A_i == coluna i  
    
    for i in range(n+m, n+m+n): # vars artificiais
        I = [0.0] * n           
        I[i-n-m] = 1.0 
        A.append(I)             # atribui por coluna e coloca mais um elem em c 
        cAUX.append(1.0)        # coloca os 1's ja prontos para passar para o formato canonico
    
    A = transpor(A)             # A_i == linha i 
    for i in range(n):          # colocar em formato canonico
        for j in range(n):
            d[i] = d[i] -1 * M[i][j] 
    for i in range(len(cAUX)):
        for j in range(n):
            cAUX[i] = cAUX[i]  -1 * A[j][i] 
    for i in range(n):
        fobj = fobj -1 * b[i]
    
    if DEBUG1: # tabela da PL aux pronta pro simplex 
        A = transpor(A) 
        printar_tabela( d, cAUX, M, A, b, fobj)    
        A = transpor(A)
        #exit(0)
    
    d, cAUX, M, A, b, fobj  = simplex(d, cAUX, M, A, b, fobj, PLaux=True) # PL aux
    
    A = transpor(A)
    A = A[:-n]      #remover as vars artificiais
    A = transpor(A)
    
    if DEBUG2: # tabela da PL aux apos o simplex  
        printar_tabela( d, c, M, A, b, fobj)
        #exit(0)
    
    return d, M, A, b
    
'''Monta o tableaux, mantem separado para facilitar o acesso.'''
def tableaux(n, m, c, A, b): # Entrada de A  [ [Restricao_i],..., [Restricao_n] ]
    d = [0.0] * n
    M = [[0.0] * n for i in range(n)] # montar matriz M de registro de opers e vetor d
    for i in range(n):
        M[i][i] = 1.0
    
    A = transpor(A) # A_i == coluna i
    dAUX = [0.0] * n
    cAUX = [0.0] * m
    
    for i in range(n):  # primeira coisa: var de folga
        I = [0.0] * n  
        I[i] = 1.0 
        A.append(I)     # atribui por coluna e coloca mais um elem em c 
        c.append(0.0)   
        cAUX.append(0.0)
    
    for i in range(n):  # segunda coisa a se fazer: b >= 0
        if b[i] < 0.0 :                 # se b_i < 0, then
            b[i] *= -1                  # flipar b_i e
            for j in range(n):          # a respectiva linha da matriz M
                M[j][i] = M[j][i] * -1
            for j in range(m+n):        # e a restricao correspondente em A
                A[j][i] = A[j][i] * -1
    
    A = transpor(A)  
    _, M, A, b = fase2(n, m, dAUX, cAUX, M, A, b) # executa PL aux e retorna uma base viavel.
    c = [i*-1 for i in c]               # flipar o vetor c, colocando o pronto para o simplex
    
    return d, c, M, A, b

'''Imprime o tableaux'''
def printar_tabela( d, c, M, A, b, fobj):
    T = []
    linha = []
    #linha = d + c + [fobj]
    linha = c + [fobj]
    T.append(linha)
    for i in range(len(b)):
        #linha = M[i] + A[i] + [b[i]]
        linha =   A[i] + [b[i]]
        T.append(linha)
    print("[", end="")
    for i in range(len(T[0])):
        if i == len(d):
            pass
            #print("|", end=" ")
        print("%.2f" % T[0][i], end="  ")
    print(" ]")
    print("[", end="")
    for i in range(len(T[0])*6):
        print("", end="-")
    print(" ]")
    for i in range(1,len(b)+1):
        print("[", end="")
        for j in range(len(T[0])):
            if j == len(d):
                pass
                #print("|", end=" ")
            print("%.2f" % T[i][j], end=", ")
        print(" ]")
    print("")

'''Imprime um vetor com 7 casas decimais'''
def printar(v):  
    for i in range(len(v)):
        f = str(Fraction(v[i]).limit_denominator())
        print("%s" % f, end=" ")
        #print("%.7f" % v[i], end=" ")
    print("")

'''Obtem o indice da menor relacao b_j/A_jk ; A_jk > 0, col e b com mesmo tamanho'''
def indice_pivo(col, b): 
    t = {}
    for i,j in enumerate(col):
        if j > 0:
            t[i] = abs(b[i])/j
    minimo = float("inf")
    indice = -1
    for i in t.keys():
        if t[i] < minimo and t[i] > 0:
            minimo = t[i]
            indice = i
    return indice

'''Executa o simplex recebendo um tableaux em forma canonica para uma base viavel e:
retorna um tableaux resolvido, caso PLaux=True
imprime o resultado caso PLaux=False'''
def simplex(d, c, M, A, b, fobj, PLaux=False): 
    transposta = False
    n = len(d)
    m = len(c)  # pega o tamanho do vetor c e consequentemente o numero de colunas de A
    cert = []
    bases = []  # lista de restricoes com bases existentes
    sol = []
    pivoLin = pivoCol = ms = -1
    
    if PLaux == False and DEBUG3:
        printar_tabela( d, c, M, A, b, fobj)    
        #exit(0)
    
    cont = 0 # este sera o contador para percorrer c, vai de 0 a m-1
    while True:
        if cont >= m: # se a condicao do final nao foi satisfeita, break
            break
        if not transposta:
            A = transpor(A) 
            transposta = True
        if c[cont] < 0:         # percorre o c procurando variaveis para otimizar
            pivoLin = indice_pivo(A[cont], b) # indice da linha do pivo
            if pivoLin == -1:   # se tiver um vetor candi
                if not transposta: 
                    A = transpor(A)
                    transposta = True
                colProblema = [-1*j for j in A[cont]] # coluna que deu problema e tornou a PL ilimitada
                if not transposta:  
                    A = transpor(A)
                    transposta = True
                for i in range(m):
                    bas = base(A[i])        # acha a solucao atraves das bases
                    if bas == -1:
                        bases.append(-1) 
                    else:
                        bases.append(bas)   # bases[i] possui o indice da coluna em A da base
                if DEBUG5:
                    if transposta:
                        A = transpor(A) 
                        transposta = False
                        printar_tabela( d, c, M, A, b, fobj) 
                        #exit(0)
                    
                for i in range(m-n):
                    if bases[i] == -1:
                        sol.append(0.0)
                        cert.append(0.0)
                    else:
                        sol.append(b[ bases[i] ])
                        cert.append(colProblema[ bases[i] ])
                    if i == cont:
                        cert[i] = 1.0 # atribuir 1 quando estiver no indice da coluna que deu problema
                print("ilimitada")
                printar(sol)
                printar(cert)
                exit(0)

            ms = A[cont][pivoLin]     # elemento do pivo, usar para zerar a coluna
            pivoCol = cont   
            fator = [1]*(1+n)         # vetor que indica por quanto temos que multiplicar cada linha para zerar a coluna do pivo.
            for i in range(n):   
                fator[i] = A[pivoCol][i] * -1   # obtem os fatores ja com o sinal correto
            fator[n] = c[pivoCol] * -1          # fator do vet c
            if transposta:   
                A = transpor(A) 
                transposta = False
            for i in range(n):      # arruma toda a linha do pivo no M e aproveita para arrumar o d
                M[pivoLin][i] = M[pivoLin][i] * (1/ms)
                d[i] = d[i] + fator[n] * M[pivoLin][i]
            
            for i in range(len(c)): # arruma toda a linha do pivo no A e aproveita para arrumar o c
                A[pivoLin][i] = A[pivoLin][i] * (1/ms)
                c[i] = c[i] + fator[n] * A[pivoLin][i]        
            b[pivoLin] = b[pivoLin] * (1/ms)
            fobj = fobj + fator[n] * b[pivoLin]
            
            for i in range(n):          # anda pelo resto das linhas de A, M e b para atualizar os valores
                if i != pivoLin :       # se a linha nao for a do pivo
                    for j in range(n):  # anda por colunas de M
                        M[i][j] = M[i][j] + fator[i] * M[pivoLin][j]
                    for j in range(len(c)): # anda pelas colunas de A
                        A[i][j] = A[i][j] + fator[i] * A[pivoLin][j]
                    b[i] = b[i] + fator[i] * b[pivoLin]
        
        cont += 1
        if cont >= m: # se andou por todos os c[cont]
            for i in range(m): 
                if c[i] < 0:    # testa se existe mais algum negativo em c
                    cont = i    # volta para o loop acima caso exista 
                    break
    
    fobj = round(fobj, 7) # aproximar numeros muito pequenos, maiores que 10e-7

    if PLaux:
        if transposta:
            A = transpor(A)
        if fobj == 0:
            return d, c, M, A, b, fobj
        elif fobj < 0:
            print("inviavel")
            printar(d)
            exit(0)
        else:
            raise ValueError("fobj PLauxiliar > 0")
    else:
        if not transposta:  
            A = transpor(A)
            transposta = True   
        for i in range(m):
            bas = base(A[i])        # obtem o indice da base
            if bas == -1:
                bases.append(-1) 
            else:
                bases.append(bas)   # bases[i] possui o indice da coluna em A da base
        
        if DEBUG4:         
            if transposta: 
                A = transpor(A) 
                printar_tabela( d, c, M, A, b, fobj)
                #exit(0)
        
        for i in range(m-n):
            if bases[i] == -1:
                sol.append(0.0)
            else:
                sol.append(b[bases[i]]) # obtem o vetor da solucao
        print("otima")
        printar([fobj])  
        print("sol")      
        printar(sol)
        printar(d)
        exit(0)

'''Resolve um problema de programacao linear no formato max cx; s.a. Ax <= b ;x >= 0'''
def main():
    convert = lambda x : list(map(float, x))
    [n, m] = convert(input().split())
    n, m = int(n), int(m)
    c = convert(input().split())
    A = []
    for i in range(n):
        A.append(convert(input().split()))
    
    A = transpor(A)  # transpor para facil acesso a linhas ou colunas
    b = A[-1]
    del A[-1]
    A = transpor(A)
    
    '''Todas as funcoes devem receber e retornarao: A := [ [Restricao_i],..., [Restricao_n] ] '''
    '''Monta a tabela e executa o fase 2 para testar inviabilidade e achar uma base viavel inicial'''
    d, c, M, A, b = tableaux(n, m, c, A, b) 
    
    '''Executa o simplex com a matriz A e M do fase 2 e printa o resultado'''
    fobj=0
    simplex(d, c, M, A, b, fobj)
    

main()

