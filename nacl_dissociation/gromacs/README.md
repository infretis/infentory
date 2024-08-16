Remember to
```bash
export OMP_NUM_THREADS=ntomp
```
where `ntomp` matches the number specified in the wmdrun section, if not it crashes

To add 10 load paths:
```bash
for i in {2..10}
do
  cp -r 1 $i
done
```

You can increase this number if you want even more load paths, but then you also need to increase the number of interfaces in infretis.toml

