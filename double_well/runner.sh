bash clean.sh
cp -r load_copy load
infretisrun -i infretis.toml
python3 change.py
bash wham.sh
