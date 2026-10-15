#include "tableau.h"

int main() {
    Tableau tab(2);
    tab.apply_h(0);
    tab.print();
    tab.apply_cx(0, 1);
    tab.print();
    tab.measure_z(1);
    tab.print();
    tab.measure_x(0);
    tab.print();
}