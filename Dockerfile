FROM python:3.11
ARG HTTPS_PROXY

WORKDIR /usr/src/app
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["gunicorn", "-w", "1", "-b", "0.0.0.0:5001", "webapp:app"]
