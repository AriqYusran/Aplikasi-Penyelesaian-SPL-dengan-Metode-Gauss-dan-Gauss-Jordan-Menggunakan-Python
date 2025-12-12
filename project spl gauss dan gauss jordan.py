import tkinter as tk
from tkinter import messagebox


def gauss(a, b):
    n = len(a)
    m = len(a[0])

    # Bentuk matriks augmented
    aug = [a[i] + [b[i]] for i in range(n)]

    # Eliminasi maju
    for k in range(n):
        pivot_row = k
        pivot_val = aug[k][k] if k < m else 0.0

        for r in range(k + 1, n):
            if k < m and abs(aug[r][k]) > abs(pivot_val):
                pivot_val = aug[r][k]
                pivot_row = r

        if k >= m or abs(pivot_val) < 1e-12:
            continue

        # Tukar baris
        aug[k], aug[pivot_row] = aug[pivot_row], aug[k]

        # Eliminasi baris di bawahnya
        for r in range(k + 1, n):
            if abs(aug[r][k]) < 1e-12:
                continue
            faktor = aug[r][k] / aug[k][k]
            aug[r] = [aug[r][c] - faktor * aug[k][c] for c in range(m + 1)]

    # Cek konsistensi
    rank_a = 0
    rank_aug = 0
    for r in range(n):
        is_all_zero_a = all(abs(aug[r][c]) < 1e-12 for c in range(m))
        is_all_zero_aug = all(abs(aug[r][c]) < 1e-12 for c in range(m + 1))

        if not is_all_zero_a:
            rank_a += 1
        if not is_all_zero_aug:
            rank_aug += 1

        if is_all_zero_a and abs(aug[r][m]) > 1e-12:
            return "tidak_ada", "Sistem tidak konsisten (tidak ada solusi)."

    if rank_a < rank_aug:
        return "tidak_ada", "Sistem tidak konsisten (tidak ada solusi)."
    elif rank_a < m:
        return "tak_hingga", "Sistem punya banyak solusi (tak hingga)."

    # Back substitution
    x = [0.0 for _ in range(m)]
    for i in range(m - 1, -1, -1):
        if abs(aug[i][i]) < 1e-12:
            return "tak_hingga", "Sistem punya banyak solusi (tak hingga)."
        sum_ax = sum(aug[i][j] * x[j] for j in range(i + 1, m))
        x[i] = (aug[i][m] - sum_ax) / aug[i][i]

    return "unik", x


def gauss_jordan(a, b):
    n = len(a)
    m = len(a[0])

    aug = [a[i] + [b[i]] for i in range(n)]

    row = 0
    for col in range(m):
        pivot_row = None
        pivot_val = 0.0
        for r in range(row, n):
            if abs(aug[r][col]) > abs(pivot_val):
                pivot_val = aug[r][col]
                pivot_row = r

        if pivot_row is None or abs(pivot_val) < 1e-12:
            continue

        # Tukar baris
        aug[row], aug[pivot_row] = aug[pivot_row], aug[row]

        # Normalkan pivot jadi 1
        pivot = aug[row][col]
        aug[row] = [val / pivot for val in aug[row]]

        # Nolkan kolom pivot pada baris lain
        for r in range(n):
            if r != row:
                faktor = aug[r][col]
                aug[r] = [aug[r][c] - faktor * aug[row][c] for c in range(m + 1)]

        row += 1
        if row == n:
            break

    # Cek konsistensi
    rank_a = 0
    rank_aug = 0
    for r in range(n):
        is_all_zero_a = all(abs(aug[r][c]) < 1e-12 for c in range(m))
        is_all_zero_aug = all(abs(aug[r][c]) < 1e-12 for c in range(m + 1))

        if not is_all_zero_a:
            rank_a += 1
        if not is_all_zero_aug:
            rank_aug += 1

        if is_all_zero_a and abs(aug[r][m]) > 1e-12:
            return "tidak_ada", "Sistem tidak konsisten (tidak ada solusi)."

    if rank_a < rank_aug:
        return "tidak_ada", "Sistem tidak konsisten (tidak ada solusi)."
    elif rank_a < m:
        return "tak_hingga", "Sistem punya banyak solusi (tak hingga)."

    solusi = [aug[i][m] for i in range(m)]
    return "unik", solusi


class SPLApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Aplikasi SPL Metode Gauss & Gauss-Jordan")

        # Frame atas: input n dan m
        top_frame = tk.Frame(root)
        top_frame.pack(pady=10)

        tk.Label(top_frame, text="Jumlah persamaan (n):").grid(row=0, column=0, padx=5, pady=5)
        self.entry_n = tk.Entry(top_frame, width=5)
        self.entry_n.grid(row=0, column=1, padx=5, pady=5)

        tk.Label(top_frame, text="Jumlah variabel (m):").grid(row=0, column=2, padx=5, pady=5)
        self.entry_m = tk.Entry(top_frame, width=5)
        self.entry_m.grid(row=0, column=3, padx=5, pady=5)

        self.btn_set = tk.Button(top_frame, text="Set Matriks", command=self.set_matrix_inputs)
        self.btn_set.grid(row=0, column=4, padx=10, pady=5)

        # Frame tengah: tempat input koefisien
        self.matrix_frame = tk.Frame(root)
        self.matrix_frame.pack(pady=10)

        # Frame pilihan metode
        method_frame = tk.Frame(root)
        method_frame.pack(pady=5)

        self.metode_var = tk.StringVar(value="gauss")
        tk.Radiobutton(method_frame, text="Metode Gauss", variable=self.metode_var, value="gauss").pack(side=tk.LEFT, padx=10)
        tk.Radiobutton(method_frame, text="Metode Gauss-Jordan", variable=self.metode_var, value="gauss_jordan").pack(side=tk.LEFT, padx=10)

        # Tombol hitung
        self.btn_hitung = tk.Button(root, text="Hitung", command=self.hitung_spl)
        self.btn_hitung.pack(pady=5)

        # Frame hasil
        result_frame = tk.Frame(root)
        result_frame.pack(pady=10)

        tk.Label(result_frame, text="Hasil:").pack(anchor="w")
        self.text_hasil = tk.Text(result_frame, width=50, height=8)
        self.text_hasil.pack()

        # Menyimpan entry matriks
        self.entries = []
        self.n = 0
        self.m = 0

    def set_matrix_inputs(self):
        # Hapus input lama
        for widget in self.matrix_frame.winfo_children():
            widget.destroy()
        self.entries = []

        try:
            self.n = int(self.entry_n.get())
            self.m = int(self.entry_m.get())
        except ValueError:
            messagebox.showerror("Error", "n dan m harus berupa bilangan bulat.")
            return

        if self.n <= 0 or self.m <= 0:
            messagebox.showerror("Error", "n dan m harus lebih dari 0.")
            return

        tk.Label(self.matrix_frame, text="Masukkan koefisien dan konstanta tiap persamaan:").grid(row=0, column=0, columnspan=self.m+1, pady=5)

        # Header
        for j in range(self.m):
            tk.Label(self.matrix_frame, text=f"a{j+1}").grid(row=1, column=j, padx=5)
        tk.Label(self.matrix_frame, text="b").grid(row=1, column=self.m, padx=5)

        # Entry matriks
        for i in range(self.n):
            row_entries = []
            for j in range(self.m + 1):  # m koefisien + 1 konstanta
                e = tk.Entry(self.matrix_frame, width=7)
                e.grid(row=2 + i, column=j, padx=3, pady=3)
                row_entries.append(e)
            self.entries.append(row_entries)

    def hitung_spl(self):
        if not self.entries:
            messagebox.showwarning("Peringatan", "Silakan set matriks terlebih dahulu.")
            return

        # Baca matriks A dan vektor b
        a = []
        b = []
        try:
            for i in range(self.n):
                row = []
                for j in range(self.m):
                    val = float(self.entries[i][j].get())
                    row.append(val)
                konst = float(self.entries[i][self.m].get())
                a.append(row)
                b.append(konst)
        except ValueError:
            messagebox.showerror("Error", "Semua koefisien dan konstanta harus berupa angka.")
            return

        metode = self.metode_var.get()
        if metode == "gauss":
            status, hasil = gauss(a, b)
        else:
            status, hasil = gauss_jordan(a, b)

        # Tampilkan hasil
        self.text_hasil.delete("1.0", tk.END)
        if status == "unik":
            for i, val in enumerate(hasil, start=1):
                self.text_hasil.insert(tk.END, f"x{i} = {val:.3f}\n")
        else:
            self.text_hasil.insert(tk.END, hasil + "\n")


if __name__ == "__main__":
    root = tk.Tk()
    app = SPLApp(root)
    root.mainloop()
