#!env /usr/bin/python3

'''

Date: 2021.11.11
Authors: hyeonjeong kim
Arranger: duaghk

Send mail module.
Refers to geninus-email.
Thanks to geninus-email authors.

'''

# import library
import smtplib
from pathlib import Path

from email import encoders
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart

from .configs import MailingConfig


class SendEmail(MailingConfig):
    """
    Auto mailing system of make_fastq results. 

    Parameters
    ----------
        Configs: dict. mailing information.
            mail_id: str, From mail id.
            mail_pw: ..
            to_mail: list. To mailing id list.

    Methods(need to arrange)
    -------
    __call__(subject, text, file, table)
        subject: str
            Mail title
        text: str
            html Mail body content with or W/O table
        file: str
            File you want to attach
        table: boolean
            A boolean flag on whether to send the file in Mail body content 

    Attribute(need to arrange)
    ---------
        from_mail: str
            Fixed    
        pw: str
            Fixed

    Example
    -------
        mail = SendEmail(to_mail='abe@abc.com,fds@dsf.com')
        mail(subject='This is a mail', text = '<html>test', file = '/data/file.excel', table = True)

    Note
    ----
    1. Only Excel file can be tabled in Mail body content 
    """

    def __init__(self, date: str):
        '''
            Input:
                mail_configs: dict. mailing information.
        '''
        MailingConfig.__init__(self)
        self.date = date

    def _connect_mail(self):
        smtp = smtplib.SMTP(host = 'gw.kr-geninus.com', port = 25)
        #smtp = smtplib.SMTP('smtp.gmail.com', 587)
        # smtp.starttls()
        smtp.login(self.login_id, self.login_pw)
        return smtp

    def make_message(self, subject, main_text, files = None):
        msg = MIMEMultipart()
        msg['Subject'] = subject
        msg['From'] = self.from_mail
        msg['To'] = ", ".join(self.to_mail_list)
        msg.attach(MIMEText(main_text, 'html'))
        if files:
            for file in files:
                attachment = open(file, 'rb')
                part = MIMEBase('application', 'octet-stream')
                part.set_payload((attachment).read())
                encoders.encode_base64(part)
                part.add_header('Content-Disposition', "attachment; filename= " + Path(file).name)
                msg.attach(part)
                # msg.attach(file)
        return msg
    
    def __call__(self, headline, main_text, filelist: list = None) -> None:
        mail_server = self._connect_mail()
        msg = self.make_message(headline, main_text, filelist)
        mail_server.send_message(msg)
        mail_server.close()
        pass










