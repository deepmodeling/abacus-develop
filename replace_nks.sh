#!/bin/bash

find ./source -type f -exec sed -i 's/\kv.nkstot_full/\kv.get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/\kv.nkstot_ibz/\kv.get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/\kv.nkstot/\kv.get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/\kv.nks\([^_[:alnum:]]\)/\kv.get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/\kv_.nkstot_full/\kv_.get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/\kv_.nkstot_ibz/\kv_.get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/\kv_.nkstot/\kv_.get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/\kv_.nks\([^_[:alnum:]]\)/\kv_.get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/kv_ptr->nkstot_full/kv_ptr->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_ptr->nkstot_ibz/kv_ptr->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_ptr->nkstot/kv_ptr->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_ptr->nks\([^_[:alnum:]]\)/kv_ptr->get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/p_kv->nkstot_full/p_kv->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/p_kv->nkstot_ibz/p_kv->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/p_kv->nkstot/p_kv->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/p_kv->nks\([^_[:alnum:]]\)/p_kv->get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/klist->nkstot_full/klist->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/klist->nkstot_ibz/klist->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/klist->nkstot/klist->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/klist->nks\([^_[:alnum:]]\)/klist->get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/pkv->nkstot_full/pkv->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/pkv->nkstot_ibz/pkv->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/pkv->nkstot/pkv->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/pkv->nks\([^_[:alnum:]]\)/pkv->get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/pkv_in->nkstot_full/pkv_in->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/pkv_in->nkstot_ibz/pkv_in->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/pkv_in->nkstot/pkv_in->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/pkv_in->nks\([^_[:alnum:]]\)/pkv_in->get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/kv->nkstot_full/kv->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv->nkstot_ibz/kv->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv->nkstot/kv->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv->nks\([^_[:alnum:]]\)/kv->get_nks\(\)\1/g' {} +

find ./source -type f -exec sed -i 's/kv_in->nkstot_full/kv_in->get_nkstot_full\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_in->nkstot_ibz/kv_in->get_nkstot_ibz\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_in->nkstot/kv_in->get_nkstot\(\)/g' {} +
find ./source -type f -exec sed -i 's/kv_in->nks\([^_[:alnum:]]\)/kv_in->get_nks\(\)\1/g' {} +
