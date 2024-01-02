FROM mysql/mysql-server:8.0

ENV TZ=UTC

#ポートを開ける
EXPOSE 3308

#MySQL設定ファイルをイメージ内にコピー
ADD ./my.cnf /etc/mysql/conf.d/my.cnf

# RUN rm -rf /var/lib/mysql

#docker runに実行される
CMD ["mysqld"]

# COPY ./mysql/my.cnf /etc/my.cnf