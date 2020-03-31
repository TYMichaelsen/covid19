# Start with Rstudio

Run the following command:

```
RSTUDIO_PASSWORD="passWORD" singularity run singularity-rstudio.simg   --auth-none 0   --auth-pam-helper rstudio_auth
```

Then login to http://<serverIP>:8787 with username of your Linux account, password as 'passWORD' or something you selected yourself.
