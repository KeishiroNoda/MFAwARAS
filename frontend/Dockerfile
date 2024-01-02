FROM node:18-alpine

WORKDIR /home/guest/SampleApp/frontend

COPY . .
RUN apk update && \
    npm install -g npm && \
    npm install -g vue-cli

CMD ["bash"]