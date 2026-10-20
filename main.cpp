#include "tableau.h"


int main() {
    Tableau tab(2);
    tab.apply_h(0);
    tab.apply_cx(0,1);
    tab.apply_z(0);
    tab.apply_cx(0, 1);
    tab.print();

    Tableau tab2(2);
    tab2.h(0);
    tab2.cx(0,1);
    tab2.z(0);
    tab2.cx(0,1);
    tab2.print();
}
