#!/bin/bash

find ./source -type f -exec sed -i 's/\.nkstot_full/\.get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/\.nkstot_ibz/\.get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/\.nkstot/\.get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/\.nks\([^_[:alnum:]]\)/\.get_nks\(\)/g' {} +

find ./source -type f -exec sed -i 's/kv_ptr->nkstot_full/kv_ptr->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_ptr->nkstot_ibz/kv_ptr->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_ptr->nkstot/kv_ptr->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_ptr->nks\([^_[:alnum:]]\)/kv_ptr->get_nks\(\)/g' {} +

find ./source -type f -exec sed -i 's/p_kv->nkstot_full/p_kv->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/p_kv->nkstot_ibz/p_kv->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/p_kv->nkstot/p_kv->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/p_kv->nks\([^_[:alnum:]]\)/p_kv->get_nks\(\)/g' {} +

find ./source -type f -exec sed -i 's/klist->nkstot_full/klist->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/klist->nkstot_ibz/klist->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/klist->nkstot/klist->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/klist->nks\([^_[:alnum:]]\)/klist->get_nks\(\)/g' {} +

