# Livro_mpc_codes
Códigos do livro Controle Preditivo Baseado em Modelo, volumes 1 e 2 para MATLAB.

**Observação 1**: ao executar os arquivos, mude a pasta de trabalho para a pasta do arquivo a ser executado.

**Observação 2**: Nos exemplos que tratam do NMPC (Volume 2, Capítulo 5), é necessário que o CasADi esteja instalado.

**Observação 3**: Algumas versões do matlab são mais sensíveis à simetria da matriz Hqp do problema de otimização. Nesses casos, logo após a linha da definição de Hqp, é possível usar o trecho de código abaixo para resolver o problema:
Hqp = (Hqp+transpose(Hqp))/2 


## Links do livro:
* [Página dos livros](https://danielml.paginas.ufsc.br/livro-mpc/)
* [Página do Volume 1 Editora Blucher](https://www.blucher.com.br/controle-preditivo-baseado-em-modelo-9788521221289)
    * [Errata Volume 1 Edição 1](https://docs.google.com/document/d/15odfF6TgR6BXdqlmBflONi6e2aLHkIWGcGHpIvZVI5A/edit?usp=sharing)
	
## Outros Links:
* [Biblioteca CasADi](https://web.casadi.org/)