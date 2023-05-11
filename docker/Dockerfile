FROM pytorch/pytorch:2.0.0-cuda11.7-cudnn8-runtime

RUN pip install \
    jupyterlab==3.6.3 \
    dask==2023.2.3

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--allow-root", "--no-browser"]
EXPOSE 8888
