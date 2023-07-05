import os
import requests

from tqdm import tqdm
from faker import Faker
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


cpu_cores = os.cpu_count()
requests.packages.urllib3.disable_warnings()


def downloader(url, type='https', show_info = True, resume=True, filename=None, num_threads=cpu_cores, timeout=10, chunk_size=1024*1000, header=None, proxies=None):
    if type == 'https':
        d = downloader_https(url, show_info, resume, filename, num_threads, timeout, chunk_size, header, proxies)
        d.download()
    else:
        pass


class downloader_https():
    def __init__(self, url, show_info = True, resume=True, filename=None, num_threads=cpu_cores, timeout=10, chunk_size=1024*1000, header=None, proxies=None):
        """
        :param url: link address
        :param filename: file name
        """
        self.url = url
        self.show_info = show_info
        self.resume = resume
        self.chunk_size = chunk_size 
        self.filename = filename
        self.num_threads = num_threads
        self.proxies = proxies
        self.timeout = timeout
        self.file_type = None
        self.accept_ranges = None
        self.content_length = None
        self.transfer_encoding = None
        if header is None:
            self.header = {}
            self.header.setdefault('User-Agent', Faker().user_agent())
        elif 'User-Agent' not in header:
            self.header.setdefault('User-Agent', Faker().user_agent())
        else:
            self.header = header

 
    def check_url(self):
        """
        check url support break-point resume and support multi-thread downloading
        """
        _, filename = os.path.split(self.url)
        self.filename = self.filename or filename

        res = requests.head(self.url, headers=self.header, proxies=self.proxies, timeout=self.timeout, allow_redirects=True, verify=False)  # verify=False 关闭ssl双向验证，解决访问https报错问题

        if not (200 <= res.status_code < 400):
            raise Exception('Bad request!')

        headers = res.headers
        self.file_type = headers.get('Content-Type')
        self.accept_ranges = headers.get('Accept-Ranges')
        self.transfer_encoding = headers.get('Transfer-Encoding')

        if self.transfer_encoding == "chunked" or self.transfer_encoding == "gzip, chunked":
            self.num_threads = 1
            self.content_length = 0
        else:
            lengths = headers.get('Content-Length')
            if lengths == None:
                self.content_length = 0
            else:
                self.content_length = int(lengths)


    def get_range(self, start=0):
        """
        set download range
        eg: [(0, 1023), (1024, 2047), (2048, 3071) ...]
        """
        if self.transfer_encoding == "chunked" or self.transfer_encoding == "gzip, chunked":
            _range = [(start, '')]
        else:
            lst = range(start, self.content_length, self.chunk_size)   
            _range = list(zip(lst[:-1], [i - 1 for i in lst[1:]]))
            _range.append((lst[-1], ''))

        return _range


    def download_by_piece(self, _range):
        start, stop = _range
        headers = {**self.header, **{"Range": f"bytes={start}-{stop}"}} # merge

        res = requests.get(self.url, headers=headers, proxies=self.proxies, timeout=self.timeout, allow_redirects=True, verify=False)
        if res.status_code != 206:
            raise Exception(f'Request raise error, url: {self.url}, range: {_range}')
        return _range, res.content


    def download(self):
        start = 0
        self.check_url()

        if self.accept_ranges != "bytes":
            if self.show_info:
                print(f'--- Mission ---: {self.url} download from scratch || with single thread, do not support breakpoint resuming')
            
            file_path = Path(self.filename)

            res = requests.get(self.url, 
                                headers=self.header, 
                                proxies=self.proxies, 
                                timeout=self.timeout, 
                                allow_redirects=True, 
                                verify=False)

            if res.status_code != 206:
                raise Exception(f'Request raise error, url: {self.url}')
 
            # write
            open(file_path, 'w').close() 
            with open(self.filename, 'rb+') as fp:
                fp.seek(0)
                fp.write(res.content)

            if self.show_info:
                print(f'--- File ---: {self.filename} download completely')
        else:
            file_path = Path(self.filename)

            if self.resume:
                open(file_path, 'w+').close()
                start = 0
                if self.show_info:
                    print(f'--- Mission ---: {self.url} download from scratch || with {self.num_threads} threads, support breakpoint resuming')

            else:
                if file_path.exists():
                    # breakpoint resume
                    start = file_path.lstat().st_size
                    if self.show_info:
                        print(f'--- Mission ---: {self.url} download from breakpoint || with {self.num_threads} threads, support breakpoint resuming')

                    # If file have already downloaded 
                    if start == self.content_length:
                        if self.show_info:
                            print(f'--- File ---: {self.filename} has already been downloaded completely')
                        return
                else:
                    open(file_path, 'w+').close()
                    start = 0
                    if self.show_info:
                        print(f'--- Mission ---: {self.url} download from scratch || with {self.num_threads} threads, support breakpoint resuming')

            # init progress bar
            if self.show_info:
                pbar = tqdm(total=self.content_length,
                        initial=start,
                        unit='B',
                        unit_scale=True,
                        desc=self.filename,
                        unit_divisor=1024)
            
            # multi-thread download
            with ThreadPoolExecutor(max_workers=self.num_threads) as pool:
                res = [pool.submit(self.download_by_piece, r) for r in self.get_range(start=start)] # or use map function

                # write
                with open(self.filename, 'rb+') as fp:
                    for item in as_completed(res):
                        _range, content = item.result()
                        start, stop = _range
                        fp.seek(start)
                        fp.write(content)
                        # update progress bar
                        if self.show_info:
                            pbar.update(self.chunk_size)

            if self.show_info:
                pbar.close()
                print(f'--- File ---: {self.filename} download completely')


    def print(self):
        self.check_url()
        print( "self.url = ", self.url, '\n',
            "self.resume = ", self.resume, '\n',
            "self.chunk_size = ", self.chunk_size, '\n',
            "self.filename = ", self.filename, '\n',
            "self.num_threads = ", self.num_threads, '\n',
            "self.proxies = ", self.proxies, '\n',
            "self.timeout = ", self.timeout, '\n',
            "self.file_type = ", self.file_type, '\n',
            "self.accept_ranges = ", self.accept_ranges, '\n',
            "self.content_length = ", self.content_length, '\n',
            "self.transfer_encoding = ", self.transfer_encoding, '\n',
            "self.header = ", self.header
            )

