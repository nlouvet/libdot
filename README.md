# libdot

Libdot library

# basic installation procedure :
autoreconf --install
./configure QDPATH=...
make

# Pour refabriquer le configure avec conviction :
make distclean
autoreconf --force --install
Comme ça, autoreconf considère que tous les fichiers sont nouveaux.

Ensuite, par exemple :

./configure CXX=g++ --enable-mpfr QDPATH=/home/nlouvet/prognlouvet/qd/qd-2.3.20-def
./configure CXX=g++-8 ARCH=avx2 QDPATH=/home/nlouvet/qd-2.3.22/ --enable-mpfr

# Dans le configure.ac:

On peut utiliser la macro AC_TRY_RUN pour lancer des tests :

AC_TRY_RUN([
#include <stdio.h>

int main(void) {
  printf("Hello!\n");
  return(0);
}
], [AC_MSG_NOTICE([cool])], [AC_MSG_FAILURE([not cool])])

Mais je ne sais pas trop comment contrôler la commande de compilation qui est utilisée, donc c'est un peut génant...
