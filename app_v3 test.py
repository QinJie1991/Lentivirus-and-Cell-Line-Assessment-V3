"""
慢病毒包装与细胞系构建评估系统 v3.4
生产级安全加固版：bcrypt密码哈希、会话超时、WAL模式、智能重试
"""

import streamlit as st
import requests
import json
import time
import re
import html
import http.client
import uuid
import logging
import sqlite3
import os
import hashlib
import bcrypt
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime, timedelta
from io import BytesIO
import pandas as pd
import urllib.request

# ==================== 安全配置 ====================
class SecurityConfig:
    """安全常量与配置"""
    GENE_NAME_PATTERN = re.compile(r'^[a-zA-Z][a-zA-Z0-9]*(-?[a-zA-Z0-9]+)*$')
    MAX_GENE_LENGTH = 50
    SESSION_TIMEOUT_MINUTES = 30  # 会话超时时间
    
    # 载体容量阈值
    STANDARD_CAPACITY = 2000
    MAX_PHYSICAL_CAPACITY = 8000
    
    # 缓存配置
    CACHE_EXPIRY_DAYS = 90

# ==================== 日志脱敏配置 ====================
class SensitiveFilter(logging.Filter):
    """过滤日志中的敏感信息"""
    def filter(self, record):
        if isinstance(record.msg, str):
            record.msg = re.sub(r'api_key=[a-zA-Z0-9]+', 'api_key=***', record.msg)
            record.msg = re.sub(r'api-key:[a-zA-Z0-9]+', 'api-key:***', record.msg, flags=re.IGNORECASE)
            record.msg = re.sub(r'Bearer [a-zA-Z0-9]+', 'Bearer ***', record.msg)
            record.msg = re.sub(r'password=.*?(?=&|\s|$)', 'password=***', record.msg)
            record.msg = re.sub(r'email=[^&\s]+', 'email=***', record.msg)
        return True

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("lentivirus_assessment")
logger.addFilter(SensitiveFilter())

# ==================== 页面配置 ====================
st.set_page_config(
    page_title="慢病毒包装与细胞系构建评估系统 v3.4",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==================== 核心数据库 ====================
class CoreDatabases:
    """核心数据库 - 高置信度基因列表"""
    
    CORE_ESSENTIAL = {
        'ACTB': ('PMID:30971823', 'DepMap核心必需', '细胞骨架结构蛋白'),
        'GAPDH': ('PMID:30971823', 'DepMap核心必需', '糖酵解关键酶'),
        'HSP90AA1': ('PMID:30971823', 'DepMap核心必需', '分子伴侣'),
        'RPL11': ('PMID:30971823', 'DepMap核心必需', '核糖体大亚基蛋白'),
        'RPS3': ('PMID:30971823', 'DepMap核心必需', '核糖体小亚基蛋白'),
        'PCNA': ('PMID:30971823', 'DepMap核心必需', 'DNA复制辅助蛋白'),
        'TOP2A': ('PMID:30971823', 'DepMap核心必需', 'DNA拓扑异构酶II'),
        'AURKB': ('PMID:30971823', 'DepMap核心必需', '有丝分裂激酶'),
        'PLK1': ('PMID:30971823', 'DepMap核心必需', '细胞周期调控激酶'),
        'BUB1': ('PMID:30971823', 'DepMap核心必需', '纺锤体检查点'),
        'CDC20': ('PMID:30971823', 'DepMap核心必需', '细胞周期后期促进复合物'),
        'CHEK1': ('PMID:30971823', 'DepMap核心必需', 'DNA损伤检查点'),
        'KIF11': ('PMID:30971823', 'DepMap核心必需', '有丝分裂驱动蛋白'),
        'PSMD1': ('PMID:30971823', 'DepMap核心必需', '蛋白酶体亚基'),
        'POLR2A': ('PMID:30971823', 'DepMap核心必需', 'RNA聚合酶II最大亚基'),
    }
    
    CORE_TOXIC = {
        'BAX': ('PMID:10625696', '促凋亡Bcl-2家族', '过表达直接激活线粒体凋亡途径'),
        'BAK1': ('PMID:10625696', '促凋亡Bcl-2家族', '线粒体外膜通透化诱导凋亡'),
        'BID': ('PMID:10625696', '促凋亡BH3-only蛋白', '连接死亡受体与线粒体凋亡'),
        'PUMA': ('PMID:12968034', 'p53下游促凋亡', '强力促凋亡BH3-only蛋白'),
        'NOXA': ('PMID:12968034', 'p53下游促凋亡', '促凋亡BH3-only蛋白'),
        'CASP3': ('PMID:9228057', '凋亡执行caspase', '过表达直接激活凋亡级联反应'),
        'CASP7': ('PMID:9228057', '凋亡执行caspase', '细胞凋亡执行分子'),
        'CASP8': ('PMID:9228057', '凋亡启动caspase', '死亡受体通路启动分子'),
        'CASP9': ('PMID:9228057', '凋亡启动caspase', '线粒体通路启动分子'),
        'FAS': ('PMID:8666142', '死亡受体', '激活外源性凋亡途径'),
        'TNF': ('PMID:15157675', '促炎细胞因子', '诱导细胞坏死性凋亡'),
        'TRAIL': ('PMID:10578115', 'TNF家族凋亡诱导配体', '选择性诱导肿瘤细胞凋亡'),
        'TP53': ('PMID:20154749', '肿瘤抑制因子', '过表达诱导G1阻滞和凋亡'),
        'CDKN1A': ('PMID:8242752', '细胞周期抑制因子', 'p21强抑制剂导致细胞周期停滞'),
        'PARP1': ('PMID:16794554', 'DNA修复酶', '过度激活导致NAD+耗竭和细胞死亡'),
    }
    
    CORE_ANTIVIRAL = {
        'MX1': ('PMID:21694717', 'ISG-I型干扰素诱导蛋白', '抑制流感病毒等RNA病毒复制'),
        'MX2': ('PMID:21694717', 'ISG-I型干扰素诱导蛋白', '抑制HIV-1等逆转录病毒'),
        'OAS1': ('PMID:21694717', 'ISG-寡腺苷酸合成酶', '激活RNase L降解病毒RNA'),
        'OAS2': ('PMID:21694717', 'ISG-寡腺苷酸合成酶', '抗病毒先天免疫效应分子'),
        'OAS3': ('PMID:21694717', 'ISG-寡腺苷酸合成酶', '抗病毒先天免疫效应分子'),
        'RNASEL': ('PMID:21694717', 'ISG-核糖核酸酶L', 'OAS通路下游降解病毒RNA'),
        'ISG15': ('PMID:21694717', 'ISG-干扰素刺激基因15', 'ISG化修饰抑制病毒复制'),
        'IFIT1': ('PMID:21694717', 'ISG-干扰素诱导蛋白', '抑制病毒翻译起始'),
        'IFIT2': ('PMID:21694717', 'ISG-干扰素诱导蛋白', '抑制病毒蛋白合成'),
        'IFIH1': ('PMID:21694717', 'ISG-MDA5', '识别病毒dsRNA激活免疫反应'),
        'DDX58': ('PMID:21694717', 'ISG-RIG-I', '识别病毒RNA诱导I型干扰素'),
        'TRIM5': ('PMID:15890885', '限制因子', '限制HIV-1等逆转录病毒复制'),
        'APOBEC3G': ('PMID:12134021', '限制因子', '胞嘧啶脱氨酶抑制HIV-1（Vif靶向）'),
        'BST2': ('PMID:19543227', '限制因子', 'Tetherin限制病毒出芽'),
    }
    
    @classmethod
    def check_gene(cls, gene_name: str, check_type: str) -> Optional[Tuple[str, str, str]]:
        gene_upper = gene_name.upper()
        if check_type == 'essential' and gene_upper in cls.CORE_ESSENTIAL:
            return cls.CORE_ESSENTIAL[gene_upper]
        elif check_type == 'toxic' and gene_upper in cls.CORE_TOXIC:
            return cls.CORE_TOXIC[gene_upper]
        elif check_type == 'antiviral' and gene_upper in cls.CORE_ANTIVIRAL:
            return cls.CORE_ANTIVIRAL[gene_upper]
        return None

# ==================== 输入验证 ====================
class InputValidator:
    @staticmethod
    def sanitize_input(text: str, max_length: int = 100) -> str:
        if not text:
            return ""
        text = text.strip()[:max_length]
        text = ''.join(char for char in text if ord(char) >= 32 and char not in ['<', '>', '"', "'"])
        return text
    
    @staticmethod
    def validate_gene_name(gene_name: str) -> Tuple[bool, str]:
        if not gene_name:
            return False, "基因名不能为空"
        if len(gene_name) > SecurityConfig.MAX_GENE_LENGTH:
            return False, f"基因名过长（最大{SecurityConfig.MAX_GENE_LENGTH}字符）"
        if not SecurityConfig.GENE_NAME_PATTERN.match(gene_name):
            return False, "基因名格式无效（必须以字母开头）"
        return True, ""
    
    @staticmethod
    def validate_password_strength(password: str) -> Tuple[bool, str]:
        if len(password) < 8:
            return False, "密码长度至少8位"
        
        has_upper = bool(re.search(r'[A-Z]', password))
        has_lower = bool(re.search(r'[a-z]', password))
        has_digit = bool(re.search(r'\d', password))
        has_special = bool(re.search(r'[!@#$%^&*(),.?":{}|<>]', password))
        
        complexity = sum([has_upper, has_lower, has_digit, has_special])
        if complexity < 3:
            return False, "密码需包含大小写字母、数字、特殊字符中的至少3种"
        
        common_weak = ['password', '123456', 'qwerty', 'admin', 'letmein']
        if password.lower() in common_weak:
            return False, "密码过于简单，请避免使用常见弱密码"
        
        return True, ""

# ==================== 增强密码管理（bcrypt哈希） ====================
class SecureAuthManager:
    """安全的认证管理器 - 使用bcrypt哈希"""
    
    @staticmethod
    def hash_password(password: str) -> str:
        """生成密码哈希"""
        salt = bcrypt.gensalt(rounds=12)
        return bcrypt.hashpw(password.encode('utf-8'), salt).decode('utf-8')
    
    @staticmethod
    def verify_password(password: str, hashed: str) -> bool:
        """验证密码"""
        try:
            return bcrypt.checkpw(password.encode('utf-8'), hashed.encode('utf-8'))
        except Exception:
            return False
    
    @staticmethod
    def check_session_timeout():
        """检查会话是否超时"""
        if "last_activity" in st.session_state:
            last_activity = st.session_state["last_activity"]
            if datetime.now() - last_activity > timedelta(minutes=SecurityConfig.SESSION_TIMEOUT_MINUTES):
                # 会话超时，清除认证状态
                for key in ["password_correct", "last_activity", "user_role"]:
                    st.session_state.pop(key, None)
                return True
        return False
    
    @staticmethod
    def update_activity():
        """更新最后活动时间"""
        st.session_state["last_activity"] = datetime.now()
    
    @staticmethod
    def logout():
        """登出功能"""
        keys_to_clear = ["password_correct", "last_activity", "user_role", "ncbi_email_input", "ncbi_key_input"]
        for key in keys_to_clear:
            st.session_state.pop(key, None)
        st.rerun()
    
    @classmethod
    def check_password(cls):
        """完整的密码验证流程"""
        # 检查会话超时
        if cls.check_session_timeout():
            st.warning("⏱️ 会话已超时（30分钟无活动），请重新登录")
        
        def password_entered():
            entered = st.session_state.get("password_input", "")
            stored_hash = st.secrets.get("APP_PASSWORD_HASH", "")
            
            # 如果没有设置密码哈希，检查明文（向后兼容）或拒绝
            if not stored_hash:
                # 尝试读取明文并自动转换（首次运行）
                plain_password = st.secrets.get("APP_PASSWORD", "")
                if plain_password:
                    # 显示警告，提示管理员更新配置
                    st.session_state["auth_error"] = "系统配置需要更新：请使用APP_PASSWORD_HASH替代APP_PASSWORD"
                else:
                    st.session_state["auth_error"] = "系统未配置访问密码"
                st.session_state["password_correct"] = False
                return
            
            # 验证密码
            if cls.verify_password(entered, stored_hash):
                st.session_state["password_correct"] = True
                st.session_state["last_activity"] = datetime.now()
                st.session_state.pop("password_input", None)
                st.session_state.pop("auth_error", None)
            else:
                st.session_state["password_correct"] = False
                st.session_state["auth_error"] = "密码错误"
        
        # 检查密码哈希配置
        if not st.secrets.get("APP_PASSWORD_HASH") and not st.secrets.get("APP_PASSWORD"):
            st.error("⚠️ 系统安全配置错误：未设置访问密码")
            st.info("请在 `.streamlit/secrets.toml` 中配置 APP_PASSWORD_HASH")
            st.code("""# 生成密码哈希的Python代码：
import bcrypt
password = "YourSecurePassword123!"
hashed = bcrypt.hashpw(password.encode(), bcrypt.gensalt())
print(hashed.decode())""", language='python')
            return False
        
        # 验证密码强度（配置时）
        if "password_strength_checked" not in st.session_state:
            plain = st.secrets.get("APP_PASSWORD", "")
            if plain:  # 如果有明文密码，检查强度
                is_strong, msg = InputValidator.validate_password_strength(plain)
                if not is_strong:
                    st.error(f"⚠️ 密码强度不足 - {msg}")
                    return False
            st.session_state["password_strength_checked"] = True
        
        # 已认证状态检查超时
        if st.session_state.get("password_correct"):
            cls.update_activity()
            return True
        
        # 显示登录界面
        st.text_input("请输入访问密码", type="password", on_change=password_entered, key="password_input")
        
        if st.session_state.get("auth_error"):
            st.error(f"😕 {st.session_state['auth_error']}")
        
        return False

# ==================== HPA数据管理（WAL模式+进度条） ====================
class HPADataManager:
    HPA_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
    LOCAL_DIR = "hpa_data"
    DB_FILE = "hpa_cache.db"
    
    def __init__(self):
        self.local_dir = self.LOCAL_DIR
        self.db_path = os.path.join(self.local_dir, self.DB_FILE)
        self.data_file = os.path.join(self.local_dir, "proteinatlas.tsv")
        self._ensure_directory()
        self._init_database()
    
    def _ensure_directory(self):
        """确保数据目录存在并设置权限"""
        if not os.path.exists(self.local_dir):
            os.makedirs(self.local_dir, mode=0o700)  # 仅所有者可读写执行
        
        # 确保目录权限正确（Unix/Linux）
        try:
            os.chmod(self.local_dir, 0o700)
        except Exception:
            pass  # Windows可能不支持
    
    def _init_database(self):
        """初始化数据库，启用WAL模式"""
        conn = sqlite3.connect(self.db_path, timeout=30.0)
        try:
            cursor = conn.cursor()
            # 启用WAL模式提高并发性能
            cursor.execute("PRAGMA journal_mode=WAL")
            cursor.execute("PRAGMA synchronous=NORMAL")
            
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS cell_line_expression (
                    gene_name TEXT,
                    cell_line TEXT,
                    rna_level TEXT,
                    protein_level TEXT,
                    reliability TEXT,
                    last_updated TIMESTAMP,
                    PRIMARY KEY (gene_name, cell_line)
                )
            ''')
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT,
                    updated_at TIMESTAMP
                )
            ''')
            # 创建索引加速查询
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_gene_cell ON cell_line_expression(gene_name, cell_line)')
            conn.commit()
        finally:
            conn.close()
    
    def check_and_download(self):
        needs_download = False
        if not os.path.exists(self.data_file):
            needs_download = True
        else:
            conn = sqlite3.connect(self.db_path, timeout=10.0)
            try:
                cursor = conn.cursor()
                cursor.execute("SELECT value, updated_at FROM metadata WHERE key='last_check'")
                result = cursor.fetchone()
                if result:
                    last_check = datetime.fromisoformat(result[1])
                    if datetime.now() - last_check > timedelta(days=30):
                        needs_download = True
                else:
                    needs_download = True
            finally:
                conn.close()
        
        if needs_download:
            self._download_hpa_data_with_progress()
    
    def _download_hpa_data_with_progress(self):
        """带进度条的下载"""
        import zipfile
        
        try:
            st.info("📥 正在下载HPA数据库（约200MB）...")
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            zip_path = os.path.join(self.local_dir, "proteinatlas.tsv.zip")
            
            # 获取文件大小
            try:
                req = urllib.request.Request(self.HPA_URL, method='HEAD')
                with urllib.request.urlopen(req, timeout=10) as response:
                    total_size = int(response.headers.get('Content-Length', 0))
            except Exception:
                total_size = 200 * 1024 * 1024
            
            if total_size == 0:
                total_size = 200 * 1024 * 1024
            
            # 流式下载
            downloaded = 0
            chunk_size = 8192
            
            response = requests.get(self.HPA_URL, stream=True, timeout=300)
            response.raise_for_status()
            
            with open(zip_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        progress = min(downloaded / total_size, 1.0)
                        progress_bar.progress(progress, text=f"已下载: {downloaded/1024/1024:.1f} MB")
            
            # 解压
            status_text.text("正在解压数据...")
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(self.local_dir)
            
            # 设置文件权限
            try:
                os.chmod(self.data_file, 0o600)
            except Exception:
                pass
            
            # 清理缓存
            self._clean_old_cache()
            
            # 更新元数据
            conn = sqlite3.connect(self.db_path, timeout=10.0)
            try:
                cursor = conn.cursor()
                cursor.execute(
                    "INSERT OR REPLACE INTO metadata (key, value, updated_at) VALUES (?, ?, ?)",
                    ('last_check', datetime.now().isoformat(), datetime.now().isoformat())
                )
                conn.commit()
            finally:
                conn.close()
            
            os.remove(zip_path)
            progress_bar.empty()
            status_text.empty()
            st.success("✅ HPA数据下载完成")
            
        except Exception as e:
            logger.error(f"HPA download error: {e}")
            st.error(f"⚠️ HPA数据下载失败: {str(e)}")
    
    def _clean_old_cache(self):
        """清理过期缓存"""
        try:
            conn = sqlite3.connect(self.db_path, timeout=10.0)
            try:
                cursor = conn.cursor()
                expiry_date = datetime.now() - timedelta(days=SecurityConfig.CACHE_EXPIRY_DAYS)
                cursor.execute(
                    "DELETE FROM cell_line_expression WHERE last_updated < ?",
                    (expiry_date.isoformat(),)
                )
                deleted = cursor.rowcount
                conn.commit()
                if deleted > 0:
                    logger.info(f"清理了 {deleted} 条过期缓存")
                    # 压缩数据库
                    cursor.execute("VACUUM")
                    conn.commit()
            finally:
                conn.close()
        except Exception as e:
            logger.error(f"Cache cleanup error: {e}")
    
    def get_expression_data(self, gene_name: str, cell_line: str) -> Optional[Dict]:
        try:
            conn = sqlite3.connect(self.db_path, timeout=10.0)
            try:
                cursor = conn.cursor()
                cursor.execute(
                    "SELECT rna_level, protein_level, reliability FROM cell_line_expression WHERE gene_name=? AND cell_line=?",
                    (gene_name.upper(), cell_line.upper())
                )
                result = cursor.fetchone()
                if result:
                    return {
                        'rna_level': result[0],
                        'protein_level': result[1],
                        'reliability': result[2],
                        'source': 'cache'
                    }
            finally:
                conn.close()
            
            if os.path.exists(self.data_file):
                return self._query_local_file(gene_name, cell_line)
            return None
        except sqlite3.Error as e:
            logger.error(f"Database error: {e}")
            return None
    
    def _query_local_file(self, gene_name: str, cell_line: str) -> Optional[Dict]:
        try:
            import csv
            gene_upper = gene_name.upper()
            cell_upper = cell_line.upper()
            cell_col_variants = [
                cell_upper.replace(" ", ""),
                cell_upper.replace(" ", "_"),
                cell_upper.replace("-", ""),
                cell_upper.replace("-", "_"),
                cell_upper
            ]
            
            with open(self.data_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                headers = reader.fieldnames
                cell_col = None
                for variant in cell_col_variants:
                    for header in headers:
                        if variant in header.upper() and ('RNA' in header.upper() or 'PROTEIN' in header.upper()):
                            cell_col = header
                            break
                    if cell_col:
                        break
                
                if not cell_col:
                    return None
                
                for row in reader:
                    if row.get('Gene', '').upper() == gene_upper or row.get('Gene name', '').upper() == gene_upper:
                        rna_key = [h for h in headers if 'RNA' in h and cell_upper.replace(" ", "") in h.upper().replace(" ", "").replace("_", "")]
                        prot_key = [h for h in headers if 'PROTEIN' in h and cell_upper.replace(" ", "") in h.upper().replace(" ", "").replace("_", "")]
                        rna_level = row.get(rna_key[0], 'Not detected') if rna_key else 'Not detected'
                        prot_level = row.get(prot_key[0], 'Not detected') if prot_key else 'Not detected'
                        result = {
                            'rna_level': rna_level,
                            'protein_level': prot_level,
                            'reliability': 'Supported',
                            'source': 'hpa_file'
                        }
                        self._cache_result(gene_name, cell_line, result)
                        return result
            return None
        except Exception as e:
            logger.error(f"HPA query error: {e}")
            return None
    
    def _cache_result(self, gene_name: str, cell_line: str, data: Dict):
        try:
            conn = sqlite3.connect(self.db_path, timeout=10.0)
            try:
                cursor = conn.cursor()
                cursor.execute('''
                    INSERT OR REPLACE INTO cell_line_expression 
                    (gene_name, cell_line, rna_level, protein_level, reliability, last_updated)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (
                    gene_name.upper(),
                    cell_line.upper(),
                    data.get('rna_level', 'Not detected'),
                    data.get('protein_level', 'Not detected'),
                    data.get('reliability', 'Unknown'),
                    datetime.now().isoformat()
                ))
                conn.commit()
            finally:
                conn.close()
        except Exception as e:
            logger.error(f"Cache write error: {e}")

# ==================== API配置（安全增强） ====================
class APIConfig:
    @staticmethod
    def get_ncbi_credentials() -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        安全获取NCBI凭证
        返回: (email, api_key, error_message)
        """
        # 优先从临时存储获取（用户输入后已清除，这里实际从session获取原始输入）
        user_email = st.session_state.get('ncbi_email_input', '').strip()
        user_key = st.session_state.get('ncbi_key_input', '').strip()
        
        # 从Secrets获取
        secret_email = st.secrets.get("NCBI_EMAIL", "")
        secret_key = st.secrets.get("NCBI_API_KEY", "")
        
        # 决定使用哪个
        email = user_email if user_email else secret_email
        api_key = user_key if user_key else secret_key
        
        # 验证
        if not email or '@' not in email:
            return None, None, "请提供有效的NCBI邮箱地址（必需）"
        
        return email, api_key, None
    
    @staticmethod
    def clear_sensitive_inputs():
        """清除敏感输入（使用后立即调用）"""
        keys = ['ncbi_key_input', 'qwen_key_input']
        for key in keys:
            if key in st.session_state:
                # 覆盖后再删除（减少内存残留时间）
                st.session_state[key] = '*' * 20
                del st.session_state[key]

# ==================== 智能重试机制（区分4xx/5xx） ====================
class RetryHandler:
    """智能重试处理器"""
    
    @staticmethod
    def retry_with_backoff(max_retries=3, backoff_factor=2, timeout=30):
        """
        装饰器：指数退避重试，智能区分客户端错误（不重试）和服务器错误（重试）
        """
        def decorator(func):
            def wrapper(*args, **kwargs):
                last_exception = None
                for attempt in range(max_retries):
                    try:
                        return func(*args, **kwargs)
                    except requests.exceptions.HTTPError as e:
                        status_code = e.response.status_code if e.response else 0
                        
                        # 4xx 客户端错误：不重试（如400参数错误、401认证失败、404未找到）
                        if 400 <= status_code < 500:
                            logger.error(f"Client error {status_code}, not retrying: {e}")
                            raise
                        
                        # 5xx 服务器错误：重试（如500内部错误、502网关错误、503服务不可用）
                        elif 500 <= status_code < 600:
                            last_exception = e
                            if attempt < max_retries - 1:
                                wait_time = backoff_factor ** attempt
                                logger.warning(f"Server error {status_code}, retrying in {wait_time}s... (attempt {attempt+1}/{max_retries})")
                                time.sleep(wait_time)
                            else:
                                raise
                        else:
                            raise
                            
                    except (requests.exceptions.Timeout, 
                            requests.exceptions.ConnectionError) as e:
                        last_exception = e
                        if attempt < max_retries - 1:
                            wait_time = backoff_factor ** attempt
                            logger.warning(f"Network error, retrying in {wait_time}s... (attempt {attempt+1}/{max_retries})")
                            time.sleep(wait_time)
                        else:
                            raise
                
                if last_exception:
                    raise last_exception
                return None
            return wrapper
        return decorator

# ==================== 频率限制器 ====================
class APIRateLimiter:
    def __init__(self, requests_per_second: float = 3.0):
        self.min_interval = 1.0 / requests_per_second
        self.last_request_time = 0
        self.lock = False  # 简单锁防止并发
    
    def wait(self):
        while self.lock:
            time.sleep(0.01)
        
        self.lock = True
        try:
            current_time = time.time()
            elapsed = current_time - self.last_request_time
            if elapsed < self.min_interval:
                time.sleep(self.min_interval - elapsed)
            self.last_request_time = time.time()
        finally:
            self.lock = False

ncbi_limiter = APIRateLimiter(3.0)

# ==================== NCBI客户端（安全增强） ====================
class NCBIClient:
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        self.email = email
        self.api_key = api_key
    
    def _make_request(self, endpoint: str, params: Dict, retmode: str = "json") -> Optional[Dict]:
        ncbi_limiter.wait()
        
        # 构建请求参数（不记录敏感信息）
        safe_params = params.copy()
        safe_params.update({
            'tool': 'LentivirusAssessment_v3',
            'email': self.email
        })
        if self.api_key:
            safe_params['api_key'] = self.api_key
        
        url = f"{self.BASE_URL}/{endpoint}"
        
        @RetryHandler.retry_with_backoff(max_retries=3, backoff_factor=2)
        def do_request():
            response = requests.get(url, params=safe_params, timeout=30)
            response.raise_for_status()
            return response.json() if retmode == "json" else response.text
        
        try:
            return do_request()
        except requests.exceptions.HTTPError as e:
            status_code = e.response.status_code if e.response else 0
            if status_code == 401 or status_code == 403:
                logger.error(f"NCBI authentication failed: {status_code}")
                raise Exception(f"NCBI API认证失败（状态码: {status_code}），请检查API Key")
            elif status_code == 429:
                logger.error("NCBI rate limit exceeded")
                raise Exception("NCBI API请求频率超限，请稍后再试或添加API Key")
            else:
                logger.error(f"NCBI request failed: {e}")
                return None
        except Exception as e:
            logger.error(f"NCBI request failed after retries: {e}")
            return None
    
    @st.cache_data(ttl=3600, show_spinner=False)
    def fetch_gene_data(_self, gene_name: str, organism: str) -> Tuple[Dict, List[Dict]]:
        search_params = {
            'db': 'gene',
            'term': f"{gene_name}[Gene] AND {organism}[Organism]",
            'retmode': 'json',
            'retmax': 1
        }
        result = _self._make_request('esearch.fcgi', search_params)
        if not result:
            return {}, []
        
        gene_ids = result.get('esearchresult', {}).get('idlist', [])
        if not gene_ids:
            return {}, []
        
        gene_id = gene_ids[0]
        summary_params = {'db': 'gene', 'id': gene_id, 'retmode': 'json'}
        result = _self._make_request('esummary.fcgi', summary_params)
        
        if not result:
            return {}, {}
        
        summary = result.get('result', {}).get(gene_id, {})
        gene_info = {
            'id': gene_id,
            'name': gene_name,
            'description': summary.get('description', ''),
            'organism': organism,
            'summary': summary.get('summary', '')
        }
        transcripts = _self._fetch_transcripts(gene_id)
        return gene_info, transcripts
    
    def _fetch_transcripts(self, gene_id: str) -> List[Dict]:
        try:
            search_params = {
                'db': 'nuccore',
                'term': f"{gene_id}[GeneID] AND (NM_[Title] OR XM_[Title])",
                'retmode': 'json',
                'retmax': 10
            }
            result = self._make_request('esearch.fcgi', search_params)
            if not result:
                return []
            
            ids = result.get('esearchresult', {}).get('idlist', [])
            if not ids:
                return []
            
            summary_params = {'db': 'nuccore', 'id': ','.join(ids), 'retmode': 'json'}
            result = self._make_request('esummary.fcgi', summary_params)
            if not result:
                return []
            
            docs = result.get('result', {})
            transcripts = []
            for uid in ids:
                try:
                    doc = docs.get(uid, {})
                    acc = doc.get('accessionversion', '')
                    length = doc.get('slen', 0)
                    if (acc.startswith('NM_') or acc.startswith('XM_')) and length > 0:
                        transcripts.append({
                            'id': acc,
                            'length': int(length),
                            'title': str(doc.get('title', ''))[:200]
                        })
                except Exception:
                    continue
            return transcripts
        except Exception as e:
            logger.error(f"Transcript fetch error: {e}")
            return []
    
    def search_gene_property_literature(self, gene_name: str, property_type: str) -> List[Dict]:
        query_map = {
            'essential': [f"{gene_name} knockout lethal", f"{gene_name} essential gene"],
            'toxic': [f"{gene_name} overexpression cytotoxic", f"{gene_name} apoptosis"],
            'antiviral': [f"{gene_name} interferon antiviral", f"{gene_name} ISG"]
        }
        
        queries = query_map.get(property_type, [f"{gene_name}"])
        all_papers = []
        seen_pmids = set()
        
        for query in queries:
            try:
                search_params = {
                    'db': 'pubmed',
                    'term': query,
                    'retmode': 'json',
                    'retmax': 5,
                    'sort': 'relevance'
                }
                result = self._make_request('esearch.fcgi', search_params)
                if not result:
                    continue
                
                pmids = result.get('esearchresult', {}).get('idlist', [])
                new_pmids = [p for p in pmids if p not in seen_pmids]
                
                if not new_pmids:
                    continue
                
                fetch_params = {'db': 'pubmed', 'id': ','.join(new_pmids), 'retmode': 'json'}
                result = self._make_request('esummary.fcgi', fetch_params)
                if not result:
                    continue
                
                docs = result.get('result', {})
                for pmid in new_pmids:
                    try:
                        doc = docs.get(pmid, {})
                        title = doc.get('title', '')
                        abstract = doc.get('abstract', '') or doc.get('sorttitle', '')
                        if not title:
                            continue
                        
                        all_papers.append({
                            'pmid': str(pmid),
                            'title': html.escape(str(title)[:300]),
                            'abstract': html.escape(str(abstract)[:800]) if abstract else "[无摘要]",
                            'query': query,
                            'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                        })
                        seen_pmids.add(pmid)
                    except Exception:
                        continue
            except Exception as e:
                logger.error(f"Literature search error: {e}")
                continue
        
        return all_papers

# ==================== 数据模型 ====================
@dataclass
class HardRuleCheck:
    rule_name: str
    passed: bool
    reason: str
    source: str
    pmid: Optional[str] = None
    pmid_list: List[str] = field(default_factory=list)
    overrideable: bool = False
    evidence_papers: List[Dict] = field(default_factory=list)
    check_level: str = "core"

# ==================== 混合硬性规则引擎 ====================
class HybridHardRulesEngine:
    def __init__(self, ncbi_client: NCBIClient):
        self.ncbi = ncbi_client
    
    def check_all(self, gene_name: str, transcripts: List[Dict], 
                  experiment_type: str) -> Tuple[bool, List[HardRuleCheck], Dict]:
        checks = []
        evidence_summary = {
            'essential_checked': False,
            'toxic_checked': False,
            'antiviral_checked': False,
            'core_hits': [],
            'literature_hits': []
        }
        
        # 载体容量检查
        if experiment_type.lower() == 'overexpression':
            check = self._check_vector_capacity(gene_name, transcripts)
            checks.append(check)
        
        # 必需基因检查
        if experiment_type.lower() == 'knockout':
            core_result = CoreDatabases.check_gene(gene_name, 'essential')
            if core_result:
                pmid, source, desc = core_result
                check = HardRuleCheck(
                    rule_name="必需基因检查（核心数据库）",
                    passed=False,
                    reason=f"❌ {gene_name}是{source}（{desc}），敲除可能导致细胞死亡",
                    source=f"DepMap数据库（{source}）",
                    pmid=pmid,
                    overrideable=False,
                    check_level="core"
                )
                checks.append(check)
                evidence_summary['essential_checked'] = True
                evidence_summary['core_hits'].append('essential')
            else:
                lit_check = self._check_by_literature(gene_name, 'essential')
                checks.append(lit_check)
                if not lit_check.passed:
                    evidence_summary['essential_checked'] = True
                    evidence_summary['literature_hits'].append('essential')
        
        # 毒性/抗病毒检查
        if experiment_type.lower() == 'overexpression':
            core_toxic = CoreDatabases.check_gene(gene_name, 'toxic')
            if core_toxic:
                pmid, source, desc = core_toxic
                check = HardRuleCheck(
                    rule_name="毒性基因检查（核心数据库）",
                    passed=False,
                    reason=f"❌ {gene_name}是{source}（{desc}），过表达可能导致细胞死亡",
                    source=f"毒性基因数据库（{source}）",
                    pmid=pmid,
                    overrideable=False,
                    check_level="core"
                )
                checks.append(check)
                evidence_summary['toxic_checked'] = True
                evidence_summary['core_hits'].append('toxic')
            else:
                lit_check = self._check_by_literature(gene_name, 'toxic')
                checks.append(lit_check)
                if not lit_check.passed:
                    evidence_summary['toxic_checked'] = True
                    evidence_summary['literature_hits'].append('toxic')
            
            core_antiviral = CoreDatabases.check_gene(gene_name, 'antiviral')
            if core_antiviral:
                pmid, source, desc = core_antiviral
                check = HardRuleCheck(
                    rule_name="抗病毒基因检查（核心数据库）",
                    passed=False,
                    reason=f"❌ {gene_name}是{source}（{desc}），过表达可能抑制慢病毒包装",
                    source=f"ISG数据库（{source}）",
                    pmid=pmid,
                    overrideable=False,
                    check_level="core"
                )
                checks.append(check)
                evidence_summary['antiviral_checked'] = True
                evidence_summary['core_hits'].append('antiviral')
            else:
                lit_check = self._check_by_literature(gene_name, 'antiviral')
                checks.append(lit_check)
                if not lit_check.passed:
                    evidence_summary['antiviral_checked'] = True
                    evidence_summary['literature_hits'].append('antiviral')
        
        return all(c.passed for c in checks), checks, evidence_summary
    
    def _check_vector_capacity(self, gene_name: str, transcripts: List[Dict]) -> HardRuleCheck:
        valid_lengths = [t.get('length', 0) for t in transcripts if t.get('length', 0) > 0]
        
        if not valid_lengths:
            return HardRuleCheck(
                rule_name="载体容量检查（过表达）",
                passed=True,
                reason="转录本长度信息暂不可获得",
                source="NCBI nuccore数据库",
                overrideable=True,
                check_level="core"
            )
        
        max_length = max(valid_lengths)
        standard = SecurityConfig.STANDARD_CAPACITY
        physical = SecurityConfig.MAX_PHYSICAL_CAPACITY
        
        if max_length <= standard:
            return HardRuleCheck(
                rule_name="载体容量检查（过表达）",
                passed=True,
                reason=f"✅ 转录本长度 {max_length}bp ≤{standard}bp，适合标准过表达",
                source="NCBI nuccore数据库",
                overrideable=True,
                check_level="core"
            )
        elif max_length <= physical:
            return HardRuleCheck(
                rule_name="载体容量检查（过表达）",
                passed=True,
                reason=f"⚠️ 转录本长度 {max_length}bp 超过{standard}bp，但低于物理极限{physical}bp，属长序列过表达",
                source="NCBI nuccore数据库",
                pmid="PMID:15819909",
                overrideable=True,
                check_level="core"
            )
        else:
            return HardRuleCheck(
                rule_name="载体容量检查（过表达）",
                passed=False,
                reason=f"❌ 转录本长度 {max_length}bp 超过物理极限{physical}bp",
                source="NCBI nuccore数据库（物理限制）",
                pmid="PMID:15819909",
                overrideable=False,
                check_level="core"
            )
    
    def _check_by_literature(self, gene_name: str, check_type: str) -> HardRuleCheck:
        papers = self.ncbi.search_gene_property_literature(gene_name, check_type)
        
        if not papers:
            type_names = {'essential': '必需性', 'toxic': '毒性', 'antiviral': '抗病毒功能'}
            return HardRuleCheck(
                rule_name=f"{type_names[check_type]}检查（文献补充）",
                passed=True,
                reason=f"✅ 核心数据库未收录，且未检索到相关文献",
                source="核心数据库+PubMed（零结果）",
                overrideable=True,
                check_level="literature"
            )
        
        evidence = []
        phrases = {
            'essential': ['lethal knockout', 'knockout is lethal', 'essential for survival'],
            'toxic': ['overexpression induced cell death', 'overexpression is cytotoxic'],
            'antiviral': ['inhibits viral replication', 'antiviral activity', 'restricts virus']
        }
        
        target_phrases = phrases.get(check_type, [])
        
        for paper in papers:
            text = (paper.get('abstract', '') + ' ' + paper.get('title', '')).lower()
            if any(phrase in text for phrase in target_phrases):
                evidence.append(paper)
        
        type_names = {'essential': '必需性', 'toxic': '毒性', 'antiviral': '抗病毒功能'}
        
        if evidence:
            pmid_list = [p['pmid'] for p in evidence[:3]]
            return HardRuleCheck(
                rule_name=f"{type_names[check_type]}检查（文献补充）",
                passed=False,
                reason=f"❌ 文献检索发现 {gene_name} 具有{type_names[check_type]}证据",
                source=f"PubMed文献检索（{len(evidence)}篇明确证据）",
                pmid=pmid_list[0],
                pmid_list=pmid_list,
                evidence_papers=evidence[:3],
                overrideable=False,
                check_level="literature"
            )
        
        return HardRuleCheck(
            rule_name=f"{type_names[check_type]}检查（文献补充）",
            passed=True,
            reason=f"✅ 检索到{len(papers)}篇文献，但未发现明确证据",
            source="PubMed文献检索（无明确证据）",
            overrideable=True,
            check_level="literature"
        )

# ==================== 基因输入组件 ====================
class GeneAutocompleteService:
    def __init__(self):
        self.clinical_tables_url = "https://clinicaltables.nlm.nih.gov/api/ncbi_genes/v3/search"
    
    @st.cache_data(ttl=3600, show_spinner=False)
    def get_suggestions(_self, query: str, organism: str = "human", limit: int = 8) -> List[Dict]:
        if not query or len(query) < 2:
            return []
        
        try:
            organism_map = {
                "human": "Homo sapiens", "mouse": "Mus musculus",
                "rat": "Rattus norvegicus", "cho": "Cricetulus griseus",
                "pig": "Sus scrofa", "monkey": "Macaca mulatta"
            }
            organism_name = organism_map.get(organism, organism)
            
            params = {
                "terms": query,
                "maxList": limit,
                "df": "symbol,name,chromosome,gene_id,type_of_gene",
                "q": f"organism:\"{organism_name}\""
            }
            
            response = requests.get(_self.clinical_tables_url, params=params, timeout=5)
            response.raise_for_status()
            data = response.json()
            
            if data and len(data) >= 3:
                results = []
                headers = data[0]
                rows = data[2]
                for row in rows:
                    gene_info = dict(zip(headers, row))
                    results.append({
                        "symbol": html.escape(gene_info.get("symbol", "")),
                        "name": html.escape(gene_info.get("name", "")),
                        "gene_id": gene_info.get("gene_id", ""),
                        "chromosome": html.escape(gene_info.get("chromosome", "")),
                        "type": html.escape(gene_info.get("type_of_gene", ""))
                    })
                return results
            return []
        except Exception as e:
            logger.warning(f"Gene suggestion error: {e}")
            return []

class GeneInputComponent:
    def __init__(self, gene_service: GeneAutocompleteService):
        self.gene_service = gene_service
    
    def render(self, organism: str, key_prefix: str = "gene") -> Optional[str]:
        input_key = f"{key_prefix}_input"
        selected_key = f"{key_prefix}_selected"
        suggestions_key = f"{key_prefix}_suggestions"
        last_query_key = f"{key_prefix}_last_query"
        
        if input_key not in st.session_state:
            st.session_state[input_key] = ""
        if selected_key not in st.session_state:
            st.session_state[selected_key] = ""
        if suggestions_key not in st.session_state:
            st.session_state[suggestions_key] = []
        if last_query_key not in st.session_state:
            st.session_state[last_query_key] = ""
        
        user_input = st.text_input(
            "基因名（支持自动完成，输入2个字符以上显示建议）",
            value=st.session_state[input_key],
            key=f"{key_prefix}_text_widget"
        )
        
        if user_input != st.session_state[input_key]:
            st.session_state[input_key] = user_input
            if st.session_state[selected_key] and user_input != st.session_state[selected_key]:
                st.session_state[selected_key] = ""
                st.session_state[suggestions_key] = []
            if len(user_input) >= 2:
                st.rerun()
        
        if len(user_input) >= 2 and not st.session_state[selected_key]:
            last_query = st.session_state.get(last_query_key, "")
            if user_input != last_query:
                suggestions = self.gene_service.get_suggestions(user_input, organism)
                st.session_state[suggestions_key] = suggestions
                st.session_state[last_query_key] = user_input
                st.rerun()
        
        suggestions = st.session_state.get(suggestions_key, [])
        if suggestions and not st.session_state[selected_key]:
            cols = st.columns(min(len(suggestions), 4))
            for i, gene in enumerate(suggestions):
                with cols[i % 4]:
                    display_text = f"{gene['symbol']}\n<small>{gene.get('name', '')[:20]}</small>"
                    if st.button(display_text, key=f"{key_prefix}_sug_{i}", use_container_width=True):
                        st.session_state[selected_key] = gene['symbol']
                        st.session_state[input_key] = gene['symbol']
                        st.session_state[f"{key_prefix}_info"] = gene
                        st.rerun()
        
        if st.session_state[selected_key] and f"{key_prefix}_info" in st.session_state:
            gene_info = st.session_state[f"{key_prefix}_info"]
            st.success(f"✅ 已选择: **{gene_info['symbol']}** | {gene_info.get('name', '')} | {gene_info.get('chromosome', '')}")
        
        if st.session_state[selected_key]:
            return st.session_state[selected_key]
        elif user_input:
            return user_input
        return None

# ==================== 报告导出 ====================
class ReportExporter:
    @staticmethod
    def generate_html_report(result: Dict) -> str:
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="UTF-8">
            <title>慢病毒包装与细胞系构建评估报告</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
                .header {{ background-color: #1f77b4; color: white; padding: 20px; text-align: center; }}
                .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
                .pass {{ color: green; }} .fail {{ color: red; }} .warning {{ color: orange; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>慢病毒包装与细胞系构建评估报告</h1>
                <p>生成时间: {result.get('timestamp', 'N/A')}</p>
            </div>
            <div class="section">
                <h2>基本信息</h2>
                <p><strong>基因:</strong> {result.get('gene', 'N/A')}</p>
                <p><strong>物种:</strong> {result.get('organism', 'N/A')}</p>
                <p><strong>细胞系:</strong> {result.get('cell_line', 'N/A')}</p>
                <p><strong>实验类型:</strong> {result.get('experiment', 'N/A')}</p>
            </div>
            <div class="section">
                <h2>评估结论</h2>
                <p>{result.get('final_recommendation', 'N/A')}</p>
                <p><strong>依据:</strong> {result.get('primary_basis', 'N/A')}</p>
            </div>
        </body>
        </html>
        """
        return html_content

# ==================== 主评估引擎 ====================
class HybridAssessmentEngine:
    def __init__(self, email: str, ncbi_api_key: Optional[str] = None):
        self.ncbi = NCBIClient(email, ncbi_api_key)
        self.hard_rules = HybridHardRulesEngine(self.ncbi)
        self.hpa = HPADataManager()
    
    def assess(self, gene_name: str, organism: str, cell_line: Optional[str], 
               experiment_type: str) -> Dict:
        result = {
            'timestamp': datetime.now().isoformat(),
            'gene': gene_name,
            'organism': organism,
            'cell_line': cell_line,
            'experiment': experiment_type,
            'decision_hierarchy': {},
            'final_recommendation': '',
            'primary_basis': ''
        }
        
        with st.spinner("🔍 检索基因详细数据..."):
            gene_info, transcripts = self.ncbi.fetch_gene_data(gene_name, organism)
        
        if not gene_info:
            return {'error': f'无法获取基因 {html.escape(gene_name)} 的信息'}
        
        result['gene_info'] = {
            'id': gene_info.get('id', ''),
            'name': gene_info.get('name', ''),
            'description': gene_info.get('description', '')[:200]
        }
        
        with st.spinner("⚙️ 执行混合硬性规则检查..."):
            hard_passed, hard_checks, evidence_summary = self.hard_rules.check_all(
                gene_name, transcripts, experiment_type
            )
        
        result['decision_hierarchy']['hard_rules'] = {
            'passed': hard_passed,
            'checks': [asdict(c) for c in hard_checks],
            'evidence_summary': evidence_summary
        }
        
        blocking = [c for c in hard_checks if not c.passed and not c.overrideable]
        if blocking:
            result['final_recommendation'] = 'BLOCKED'
            result['primary_basis'] = '硬性生物学限制（核心数据库或文献证据）'
            result['blocking_evidence'] = [asdict(c) for c in blocking]
            return result
        
        if organism == 'Homo sapiens' and cell_line:
            with st.spinner("🧬 查询HPA表达数据..."):
                hpa_data = self.hpa.get_expression_data(gene_name, cell_line)
                result['hpa_data'] = hpa_data or {
                    'rna_level': '数据难以获得',
                    'protein_level': '数据难以获得',
                    'reliability': 'N/A'
                }
        
        if cell_line:
            with st.spinner("🧫 检索细胞相关参数..."):
                cell_params = self.ncbi.search_gene_property_literature(cell_line, 'lentivirus')
                same_cell_studies = self.ncbi.search_gene_property_literature(f"{gene_name} {cell_line}", 'study')
                result['cell_assessment'] = {
                    'lentivirus_params': cell_params if cell_params else '无已报道的参数',
                    'same_cell_gene_studies': same_cell_studies if same_cell_studies else '无同细胞同基因研究报道'
                }
        
        warning_checks = [c for c in hard_checks if not c.passed and c.overrideable]
        if warning_checks:
            result['final_recommendation'] = "⚠️ 警告：检测到潜在风险，建议谨慎操作"
            result['primary_basis'] = f"基于{len(warning_checks)}项警告（可人工覆盖）"
        else:
            result['final_recommendation'] = "✅ 未检测到明确风险，可进行标准流程"
            result['primary_basis'] = "基于核心数据库筛查和文献检索"
        
        return result

# ==================== UI渲染 ====================
def render_sidebar():
    with st.sidebar:
        st.header("⚙️ API配置")
        st.subheader("NCBI配置")
        ncbi_email = st.text_input("NCBI邮箱", value="", key="ncbi_email_input", help="优先使用此处输入，留空使用Secrets")
        ncbi_key = st.text_input("NCBI API Key", type="password", key="ncbi_key_input", help="可选，用于提高访问频率限制")
        
        email, key, error = APIConfig.get_ncbi_credentials()
        if error:
            st.error(error)
        else:
            st.success("✅ NCBI API有效")
        
        st.divider()
        
        # 登出按钮
        if st.session_state.get("password_correct"):
            if st.button("🚪 登出", type="secondary"):
                SecureAuthManager.logout()
        
        st.divider()
        st.caption("🔒 核心列表+文献补充混合策略（v3.4安全加固版）")

def render_main_panel():
    st.markdown("""
    <h1 style='text-align: center; color: #1f77b4; margin-bottom: 30px;'>
        🧬 慢病毒包装与细胞系构建评估系统
    </h1>
    """, unsafe_allow_html=True)
    
    st.markdown("### 🔬 实验参数输入")
    col1, col2 = st.columns(2)
    
    with col1:
        organism = st.selectbox(
            "物种",
            ["human", "mouse", "rat", "cho", "pig", "monkey"],
            format_func=lambda x: {
                "human": "人类 (Homo sapiens)",
                "mouse": "小鼠 (Mus musculus)",
                "rat": "大鼠 (Rattus norvegicus)",
                "cho": "CHO (Cricetulus griseus)",
                "pig": "家猪 (Sus scrofa)",
                "monkey": "猴子 (Macaca mulatta)"
            }.get(x, x)
        )
        gene_service = GeneAutocompleteService()
        gene_component = GeneInputComponent(gene_service)
        gene = gene_component.render(organism, key_prefix="main_gene")
    
    with col2:
        cell_line = st.text_input("细胞名（可选）", placeholder="例如：HEK293T, HeLa, A549", 
                                 help="输入细胞系名称以获取特定细胞系的感染参数和HPA表达数据")
        exp_type = st.selectbox(
            "评估选项",
            ["overexpression", "knockdown", "knockout"],
            format_func=lambda x: {
                "overexpression": "过表达 (OE)",
                "knockdown": "敲低 (RNAi)",
                "knockout": "敲除 (CRISPR)"
            }.get(x, x)
        )
    
    analyze = st.button("🚀 开始AI智能评估", type="primary", use_container_width=True)
    return organism, gene, cell_line, exp_type, analyze

def render_results(result: Dict):
    if 'error' in result:
        st.error(result['error'])
        return
    
    st.divider()
    st.markdown(f"## 🎯 评估报告 - {html.escape(result['gene'])}")
    
    col_exp1, col_exp2 = st.columns(2)
    with col_exp1:
        if st.button("📄 导出HTML报告"):
            exporter = ReportExporter()
            html_report = exporter.generate_html_report(result)
            st.download_button(
                "下载HTML",
                html_report,
                file_name=f"assessment_{result['gene']}_{datetime.now().strftime('%Y%m%d')}.html",
                mime="text/html"
            )
    
    rec = result['final_recommendation']
    rec_color = {"❌": "#ffebee", "⚠️": "#fff3e0", "⚡": "#fff8e1", "✅": "#e8f5e9"}.get(rec[:2], "#f5f5f5")
    
    st.markdown(f"""
    <div style='padding: 20px; background-color: {rec_color}; 
                border-radius: 10px; text-align: center; margin: 20px 0;
                border: 2px solid #ddd;'>
        <h3>{html.escape(rec)}</h3>
        <small>{html.escape(result.get('primary_basis', ''))}</small>
    </div>
    """, unsafe_allow_html=True)
    
    tabs = st.tabs(["硬性规则检查", "HPA表达数据", "细胞评估", "序列设计"])
    
    with tabs[0]:
        st.markdown("### 🚦 混合硬性规则检查（核心数据库+文献补充）")
        hierarchy = result.get('decision_hierarchy', {})
        hard_rules = hierarchy.get('hard_rules', {})
        evidence_summary = hard_rules.get('evidence_summary', {})
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("核心数据库命中", len(evidence_summary.get('core_hits', [])))
        with col2:
            st.metric("文献补充发现", len(evidence_summary.get('literature_hits', [])))
        with col3:
            st.metric("总检查项", len(hard_rules.get('checks', [])))
        
        st.divider()
        
        for check in hard_rules.get('checks', []):
            icon = "✅" if check['passed'] else "❌" if not check['overrideable'] else "⚠️"
            color = "green" if check['passed'] else "red" if not check['overrideable'] else "orange"
            level_badge = "🔴 核心" if check.get('check_level') == "core" else "🔵 文献"
            
            evidence_md = ""
            if check.get('evidence_papers'):
                evidence_md = "<br/><small>📚 证据文献：</small><br/>"
                for paper in check['evidence_papers'][:2]:
                    pmid = paper.get('pmid', '')
                    title = paper.get('title', '')[:80]
                    evidence_md += f'<small>• <a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">PMID:{pmid}</a> {html.escape(title)}...</small><br/>'
            
            pmid_list = check.get('pmid_list', [])
            pmid_badge = ""
            if pmid_list:
                pmid_badge = f'<br/><small>PMID: {", ".join([f"<a href=\'https://pubmed.ncbi.nlm.nih.gov/{p}/\' target=\'_blank\'>{p}</a>" for p in pmid_list[:3]])}</small>'
            
            st.markdown(f"""
            <div style='padding: 15px; border-left: 4px solid {color}; 
                        background-color: #f8f9fa; margin: 10px 0; border-radius: 5px;'>
                <h4>{icon} {html.escape(check['rule_name'])} {level_badge}</h4>
                <p>{html.escape(check['reason'])}</p>
                <small>来源: {html.escape(check['source'])}</small>
                {pmid_badge}
                {evidence_md}
            </div>
            """, unsafe_allow_html=True)
    
    with tabs[1]:
        st.markdown("### 🧬 HPA表达量数据")
        hpa_data = result.get('hpa_data')
        if hpa_data:
            if 'rna_level' in hpa_data:
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("RNA水平", hpa_data.get('rna_level', 'N/A'))
                with col2:
                    st.metric("蛋白水平", hpa_data.get('protein_level', 'N/A'))
                with col3:
                    st.metric("可靠性", hpa_data.get('reliability', 'N/A'))
                gene = result['gene']
                st.markdown(f"[查看HPA详情](https://www.proteinatlas.org/{gene}-{result.get('cell_line', '')})")
            else:
                st.info(hpa_data.get('message', '数据难以获得'))
        else:
            st.info("⚠️ 仅当物种为人类且输入细胞名时显示HPA数据")
    
    with tabs[2]:
        st.markdown("### 🧫 细胞系构建评估数据")
        cell_data = result.get('cell_assessment')
        if cell_data:
            st.subheader("📚 同细胞同基因文献")
            studies = cell_data.get('same_cell_gene_studies')
            if studies and studies != '无同细胞同基因研究报道':
                for study in studies[:5]:
                    st.markdown(f"""
                    - **{html.escape(study.get('title', ''))}**  
                      {html.escape(study.get('journal', ''))} ({study.get('year', '')})  
                      [PMID: {study.get('pmid', '')}]({study.get('url', '')})
                    """)
            else:
                st.warning(studies if isinstance(studies, str) else "无已报道的研究")
            
            st.subheader("🦠 慢病毒MOI参数")
            params = cell_data.get('lentivirus_params')
            if params and params != '无已报道的参数':
                st.json(params)
            else:
                st.warning("无已报道的参数")
        else:
            st.info("⚠️ 输入细胞名以获取细胞评估数据")
    
    with tabs[3]:
        st.info("⚠️ 序列设计功能在基础版本中通过文献检索实现")

def main():
    # 检查认证（包含超时检查）
    if not SecureAuthManager.check_password():
        st.stop()
    
    # 更新活动时间
    SecureAuthManager.update_activity()
    
    # 初始化HPA
    hpa_manager = HPADataManager()
    hpa_manager.check_and_download()
    
    # 渲染界面
    render_sidebar()
    organism, gene, cell_line, exp_type, analyze = render_main_panel()
    
    if analyze:
        if not gene:
            st.error("⚠️ 请输入或选择一个基因")
            return
        
        is_valid, error_msg = InputValidator.validate_gene_name(gene)
        if not is_valid:
            st.error(f"输入验证失败: {error_msg}")
            return
        
        gene_clean = InputValidator.sanitize_input(gene, 50)
        cell_clean = InputValidator.sanitize_input(cell_line, 100) if cell_line else None
        
        organism_map = {
            "human": "Homo sapiens", "mouse": "Mus musculus",
            "rat": "Rattus norvegicus", "cho": "Cricetulus griseus",
            "pig": "Sus scrofa", "monkey": "Macaca mulatta"
        }
        organism_clean = organism_map.get(organism, organism)
        
        email, ncbi_key, error = APIConfig.get_ncbi_credentials()
        if error:
            st.error(error)
            return
        
        try:
            engine = HybridAssessmentEngine(email=email, ncbi_api_key=ncbi_key)
            
            with st.spinner("正在进行混合策略评估（核心数据库+文献补充）..."):
                result = engine.assess(gene_clean, organism_clean, cell_clean, exp_type)
            
            render_results(result)
            
        except Exception as e:
            error_id = str(uuid.uuid4())[:8]
            logger.exception(f"Unhandled error: {e}")
            st.error(f"❌ 系统错误（ID: {error_id}），请联系管理员")

if __name__ == "__main__":
    main()
